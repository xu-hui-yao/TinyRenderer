#include <parse/parser.h>
#include <components/scene.h>
#include <render/block.h>
#include <filesystem/resolver.h>
#include <thread>
#include <iostream>
#include <components/camera.h>
#include <core/timer.h>
#include <components/bitmap.h>
#include <vector>

using namespace tiny_renderer;

static void render_block(const std::shared_ptr<Scene>& scene, const std::shared_ptr<Sampler>& sampler,
                         ImageBlock& block) {
	const auto& camera = scene->get_camera();
	const auto& integrator = scene->get_integrator();

	Point2i offset = block.get_offset();
	Vector2i size = block.get_size();

	/* Clear the block contents */
	block.clear();

	/* For each pixel and pixel sample sample */
	for (int y = 0; y < size.y(); ++y) {
		for (int x = 0; x < size.x(); ++x) {
			for (uint32_t i = 0; i < sampler->get_sample_count(); ++i) {
				Point2f pixel_sample = Point2f(static_cast<float>(x + offset.x()), static_cast<float>(y + offset.y())) +
					sampler->next2d();
				Point2f aperture_sample = sampler->next2d();

				/* Sample a ray from the camera */
				Ray3f ray;
				bool valid = true;
				Color3f value = camera->sample_ray(ray, pixel_sample, aperture_sample, valid);

				/* Compute the incident radiance */
				value *= integrator->li(scene, sampler, ray, valid);

				/* Store in the image block */
				block.put(pixel_sample, value);
			}
		}
	}
}

static void render(const std::shared_ptr<Scene>& scene, const std::string& filename, int thread_count) {
	auto camera = scene->get_camera();
	Vector2i output_size = camera->get_output_size();
	scene->get_integrator()->preprocess(scene);

	/* Create a block generator (i.e. a work scheduler) */
	BlockGenerator block_generator(output_size, M_BLOCK_SIZE);

	/* Allocate memory for the entire output image and clear it */
	ImageBlock result(output_size, camera->get_reconstruction_filter());
	result.clear();

	/* Create threads */
	std::vector<std::thread> threads;
	int blocks_per_thread = block_generator.get_block_count() / thread_count;
	int remaining_blocks = block_generator.get_block_count() % thread_count;

	std::cout << "Rendering .. " << std::endl;
	Timer timer;

	threads.reserve(thread_count);
	for (int i = 0; i < thread_count; ++i) {
		threads.emplace_back(
			[&, i] {
				// Each thread works on its own block range
				int start_block = i * blocks_per_thread;
				int end_block = i == thread_count - 1
					                ? start_block + blocks_per_thread + remaining_blocks
					                : start_block + blocks_per_thread;

				ImageBlock block(Vector2i(M_BLOCK_SIZE), camera->get_reconstruction_filter());
				std::shared_ptr sampler(scene->get_sampler()->clone());

				for (int j = start_block; j < end_block; ++j) {
					block_generator.next(block);
					sampler->prepare(block);  // The sampler of different block has different seeds.

					// /* Render all contained pixels */
					render_block(scene, sampler, block);

					/* The image block has been processed. Now add it to
					   the "big" block that represents the entire image */
					result.put(block);
				}
			});
	}

	// Join all threads
	for (auto& t : threads) {
		t.join();
	}

	std::cout << "done. (took " << timer.elapsed_string() << ")" << std::endl;

	/* Now turn the rendered image block into
	   a properly normalized bitmap */
	auto bitmap = result.to_bitmap();

	/* Determine the filename of the output bitmap */
	filesystem::path file = get_file_resolver()->resolve(filesystem::path(filename));
	std::string output_name = file.make_absolute().str();
	size_t lastdot = output_name.find_last_of('.');

	if (lastdot != std::string::npos)
		output_name.erase(lastdot, std::string::npos);

	/* Save using the OpenEXR format */
	// bitmap->save_exr(output_name + ".exr");

	/* Save tone mapped (sRGB) output using the PNG format */
	bitmap->save_png(output_name + ".png");
}

int main(int argc, char** argv) {
	int thread_count = 1;

	if (argc < 2) {
		std::cerr << "Syntax: " << argv[0] << " <scene.xml> [--no-gui] [--threads N]" << std::endl;
		return -1;
	}

	std::string scene_name;

	for (int i = 1; i < argc; ++i) {
		std::string token(argv[i]);
		if (token == "-t" || token == "--threads") {
			if (i + 1 >= argc) {
				std::cerr << "\"--threads\" argument expects a positive integer following it." << std::endl;
				return -1;
			}
			char* end_ptr;
			thread_count = strtol(argv[i + 1], &end_ptr, 10);
			i++;
			if (thread_count <= 0) {
				std::cerr << "\"--threads\" argument expects a positive integer following it." << std::endl;
				return -1;
			}

			continue;
		}
#ifdef M_DEBUG
		thread_count = 1;
#endif

		filesystem::path path(argv[i]);

		if (path.extension() == "xml") {
			scene_name = argv[i];

			/* Add the parent directory of the scene file to the
			   file resolver. That way, the XML file can reference
			   resources (OBJ files, textures) using relative paths */
			get_file_resolver()->prepend(path.parent_path());
		}
	}

	if (scene_name.empty()) {
		std::cerr << "Please provide xml" << std::endl;
		return -1;
	}

	try {
		auto root(load_from_xml(scene_name));
		if (root->get_class_type() == Object::EScene) {
			root->construct();
			render(std::dynamic_pointer_cast<Scene>(root), scene_name, thread_count);
		}
	} catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		return -1;
	}

	return 0;
}