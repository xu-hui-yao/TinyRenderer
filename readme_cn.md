# Tiny Renderer



![Windows Build Status](https://github.com/xu-hui-yao/TinyRenderer/actions/workflows/windows.yml/badge.svg)![Linux Build Status](https://github.com/xu-hui-yao/TinyRenderer/actions/workflows/linux.yml/badge.svg)![macOS Build Status](https://github.com/xu-hui-yao/TinyRenderer/actions/workflows/macos.yml/badge.svg)



## 介绍

这是一个用 C++ 实现的路径追踪渲染器，设计遵循 [PBRT](https://pbr-book.org/3ed-2018/contents) 标准，并支持多线程渲染。项目结构参考 [Mitsuba 3](http://www.mitsuba-renderer.org/)。

该渲染器接受描述场景的 XML 文件作为输入，生成图像并保存至指定路径。

- 加速结构
  - 支持层次包围盒（BVH）、KD树（KDTree）和八叉树（Octree）用于光线求交加速。
- 材质支持
  - 基本材质：漫反射、导体、介电体、塑料、超薄介电体、粗糙导体、粗糙介电体、粗糙塑料。
  - 特殊材质：掩膜、凹凸贴图、双面材质。
  - 纹理支持：各向异性纹理（exr、hdr、png、jpeg等格式）、常量纹理、棋盘格纹理。
- 光照模型
  - 支持任意多个面光源以及天空盒光照。
- 几何形状
  - 支持 OBJ 文件定义的网格
  - 内置几何体：立方体、长方形、球体
- 路径追踪
  - 使用多重重要性采样，支持对材质以及场景中的光源进行重要性采样。
- 降噪
  - 单帧蒙特卡洛降噪器：`atrous`（保边 À-Trous 小波）、`nlm`（非局部均值）、`outlier`（萤火虫/离群点剔除）。
  - 通过场景 XML 配置，并且支持嵌套：一个降噪器可以包含一个子降噪器作为预滤波（典型用法是 `outlier` 接在 `atrous` 之前）。
  - 在色调映射之前作用于线性辐射值，因此同一套实现同时服务于 CPU 与 GPU 渲染路径。
- GPU 后端（实验性）
  - 基于 Vulkan compute 的路径追踪器，提供两种后端：**Megakernel**（每个 spp 一次 dispatch，整个弹射循环在单个 shader 内）与 **Wavefront**（每次弹射拆成 6 个 kernel，通过压缩工作队列连接）。
  - 着色器使用 [Slang](https://shader-slang.org/) 编写，构建时交叉编译为 SPIR-V。
- 实时交互模式
  - 常驻的 GPU 会话，每帧渲染 1 spp，配合 SVGF 时空滤波、轨道/飞行相机以及 ImGui 参数面板，详见[实时交互模式](#实时交互模式)。



## 编译运行

### 环境依赖

| 依赖 | 是否必需 | 用途 |
| --- | --- | --- |
| CMake >= 3.19 与支持 C++17 的编译器 | 必需 | 全部功能 |
| [Vulkan SDK](https://vulkan.lunarg.com/)（其中包含 `slangc`，且需在 `PATH` 中） | 可选 | GPU 后端（`--gpu`）与实时交互（`--interactive`） |
| [GLFW 3](https://www.glfw.org/) + OpenGL | 可选 | 预览窗口（`--progress`）与实时交互 |

第三方库（`pugixml`、`stb`、`tinyexr`、Dear ImGui、GLFW）都以 git submodule 的形式放在 `ext/` 下，无需额外安装（clone 时加 `--recursive`，或执行 `git submodule update --init --recursive`）。

注意"可选"不等于"不装也能顺利构建"：`M_ENABLE_GPU` 默认为 `ON`，所以默认构建**确实需要** Vulkan SDK 与 `slangc`；只有 GLFW 是 CMake 能在缺失时静默降级的那一个——即便如此，在 Windows/Apple 上它也会退回到编译自带的 `ext/glfw` 源码；Linux 仍然需要系统包（GLFW 的 X11/Wayland 后端需要 dev 头文件），所以只有 Linux 上"缺 GLFW"才等于"不带窗口"。

可选依赖可以用系统的包管理器安装：

```bash
# macOS
brew install glfw            # 预览窗口 / 实时交互
# 然后从 https://vulkan.lunarg.com/ 安装 LunarG Vulkan SDK（提供 Vulkan 与 slangc）

# Ubuntu / Debian
sudo apt install libglfw3-dev
# 然后安装 LunarG Vulkan SDK

# Windows：安装 Vulkan SDK（提供 Vulkan 与 slangc）；GLFW 无需安装，会直接编译自带的 ext/glfw
```

### CMake 选项

| 选项 | 默认值 | 作用 |
| --- | --- | --- |
| `M_ENABLE_GPU` | `ON` | 构建 Vulkan GPU 后端并链接进 `tiny-renderer`。需要 Vulkan SDK 与 `slangc`，**缺任意一个都会直接配置失败**（`find_package(Vulkan REQUIRED)` / `find_program(slangc REQUIRED)`）。不想要 GPU 就设 `-DM_ENABLE_GPU=OFF`。 |
| `M_ENABLE_PREVIEW_GUI` | `ON` | 构建 GLFW + Dear ImGui 窗口。如果**找不到** `glfw3`/OpenGL，CMake 只会打印一条状态信息并继续构建——此时 `--progress` 只剩下控制台进度条，而 `--interactive` 不可用。不过在 Windows/Apple 上这种情况很少发生：找不到系统 `glfw3` 时会改为编译自带的 `ext/glfw` 源码。 |

例如，一个不带窗口的纯 CPU 构建：

```bash
cmake -DCMAKE_BUILD_TYPE=Release -DM_ENABLE_GPU=OFF -DM_ENABLE_PREVIEW_GUI=OFF -S . -B build
```

### Windows 11

编译器：Visual Studio 2022

进入项目根目录，打开 cmd 运行以下指令：

```cmd
mkdir build
cd build
cmake .. -G "Visual Studio 17 2022"
cmake --build . --config Release
```

可执行文件生成于`build/src/Release`，运行可执行文件：

```cmd
tiny-renderer.exe 'xml relative path of the root directory' -t 'thread count'
```

即可在 xml 同级目录下生成渲染图（png）。

### MacOS

编译器：xcode。确保你的系统已经安装了 CMake 和 Xcode 命令行工具。

进入项目根目录，打开终端运行以下指令：

```bash
mkdir build
cd build
cmake .. -G "Xcode"
xcodebuild -configuration Release
```

可执行文件生成于`build/src/Release`，运行可执行文件：

```bash
tiny-renderer 'xml relative path of the root directory' -t 'thread count'
```

即可在 xml 同级目录下生成渲染图（png）。

### Linux

编译器：g++。

进入项目根目录，打开终端运行以下指令：

```bash
cmake -DCMAKE_BUILD_TYPE=Release -S . -B build
cd build
make -j${proc}
```

可执行文件生成于`build/src/`，运行可执行文件：

```bash
tiny-renderer 'xml relative path of the root directory' -t 'thread count'
```

即可在 xml 同级目录下生成渲染图（png）。

### 命令行参数

```
tiny-renderer <scene.xml> [options]
```

| 参数 | 说明 |
| --- | --- |
| `<scene.xml>` | 必需。场景文件路径，其所在目录会被加入文件搜索路径，因此 XML 中可以用相对路径引用 OBJ 与贴图。 |
| `-t N` / `--threads N` | CPU 渲染线程数（默认 1）。`--gpu` 与 `--interactive` 下无效。 |
| `--gpu[=megakernel\|wavefront]` | 使用 GPU 渲染，默认后端为 `megakernel`。 |
| `--progress` | 打印实时控制台进度条（百分比/已用时/预计剩余），若编译时带了预览窗口，还会弹出窗口显示当前图像。 |
| `--tonemap=none\|aces` | 输出 PNG 时使用的色调映射（默认 `none`）。 |
| `--denoise[=outlier\|atrous\|nlm]` | 启用降噪器（使用默认参数），会覆盖场景 XML 中的 `<denoiser>` 配置。只写 `--denoise` 时选择 `atrous`。 |
| `--dump-aov` | 额外把降噪器的输入缓冲导出为 PNG，便于调参。 |
| `--interactive[=WxH]` | 打开实时交互窗口。`WxH` 覆盖场景的输出分辨率（默认使用 XML 中配置的分辨率）。 |

`--interactive` 优先于 `--gpu`；交互模式**不会**输出任何 PNG，只做显示（关闭窗口即退出）。

#### 输出文件

以 `assets/teapot/teapot.xml` 为例，输出写在场景文件同级目录下：

| 文件 | 何时生成 |
| --- | --- |
| `teapot.png` | 总是生成——最终图像（若启用降噪则为降噪后的结果） |
| `teapot_noisy.png` | 只要跑了降噪器就会生成——未降噪的原图，便于 A/B 对比 |
| `teapot_albedo.png`、`teapot_normal.png`、`teapot_variance.png` | 仅在使用 `--dump-aov` 时生成 |

#### 示例

```bash
# CPU 渲染，8 线程
tiny-renderer assets/teapot/teapot.xml -t 8

# GPU megakernel，带实时预览窗口
tiny-renderer assets/dragon/dragon.xml --gpu=megakernel --progress

# GPU wavefront，降噪 + ACES 色调映射，并导出降噪器输入
tiny-renderer assets/box/box.xml --gpu=wavefront --denoise=atrous --tonemap=aces --dump-aov

# 以 1280x720 打开实时交互窗口
tiny-renderer assets/dragon/dragon.xml --interactive=1280x720
```



## 实时交互模式

`--interactive` 用**常驻的 GPU 会话**取代了"渲染一次并写出 PNG"的离线流程：每帧只 dispatch `spp_per_frame` 个样本，由 SVGF（时域重投影 + À-Trous 空域滤波）把 1 spp 的信号逐步累积成收敛图像。结果直接从 GPU 的 RGBA8 缓冲上传显示，主机端每帧不做任何逐像素处理。

### 编译前提

只有**同时**满足以下条件时，交互模式才会被编译进来（见 [CMake 选项](#cmake-选项)）：

- `M_ENABLE_GPU=ON`（Vulkan SDK + `slangc`），且
- `M_ENABLE_PREVIEW_GUI=ON`，并且在配置阶段确实找到了 `glfw3` + OpenGL（Windows/Apple 上找不到系统 GLFW 时会改用自带的 `ext/glfw`）。

否则运行 `--interactive` 会打印 `This build lacks interactive support ...` 并退出。配置时要留意 CMake 是否输出了 `M_ENABLE_PREVIEW_GUI is ON but glfw3/OpenGL were not found`——这就是"静默退化"的情况。

### 操作方式

| 输入 | 轨道模式（Orbit，默认） | 飞行模式（Fly） |
| --- | --- | --- |
| 左键拖拽 | 绕目标点旋转 | - |
| 中键/右键拖拽 | 平移 | 平移 |
| 滚轮 | 推拉镜头 | 前进/后退 |
| `W` `A` `S` `D` | - | 沿视轴移动 |
| `Q` / `E` | - | 下移 / 上移 |
| `Shift` | - | 4 倍加速 |

初始视角与 FOV 都从场景 XML 读取，因此打开窗口时看到的正是 XML 中配置的机位。相机或渲染参数一旦改变，时域累积就会重置（历史对新的视角是过期的）；仅影响显示的设置不会触发重置。

### 参数面板

左上角 ImGui 面板提供：

- **Stats**：FPS / 帧时间、已累积 spp，以及各 pass 的 GPU 耗时分解（`path_trace`、`temporal`、`atrous`、`display`）。
- **Camera**：Orbit / Fly、FOV、移动速度、重置视角。
- **Render**：每帧 spp（1-8）、最大深度、俄罗斯轮盘赌深度、辐射亮度钳制（0 表示关闭）、反照率解调、重置累积。
- **Denoise**：时域滤波开关与拒绝阈值（`alpha`、`phi_depth`、`phi_normal`、历史钳制与增长系数），以及空域 pass（迭代次数、`phi_color`、`phi_normal`、`phi_depth`）。把迭代次数设为 0 可以单独观察时域 pass 的效果。
- **Display**：曝光、色调映射（Clamp / Reinhard / ACES）、sRGB 传输函数、gamma，以及**调试视图**（raw illum、filtered illum、albedo、normal、depth、variance、历史长度、position）——调滤波参数主要靠这些视图。

### 自动化跑帧

设置环境变量 `M_INTERACTIVE_MAX_FRAMES=N` 可在跑满 N 帧后自动退出并打印平均 FPS。否则循环会一直运行到窗口被关闭，帧率也就取决于何时关闭窗口，不适合做对比。

```bash
M_INTERACTIVE_MAX_FRAMES=300 tiny-renderer assets/dragon/dragon.xml --interactive=1280x720
```



## GPU 后端

GPU 路径（`--gpu`、`--interactive`）复用与 CPU 完全相同的 `Scene`/BVH 构建流程，再把场景展平成纯数组（`Scene::build_gpu_scene()`，见 `include/core/gpu_scene.h`），一次性上传，然后用 Slang 编写、构建期交叉编译为 SPIR-V 的 Vulkan compute shader（`src/gpu/shaders`）执行路径追踪。

- **Megakernel**（`--gpu=megakernel`）：每个 spp 一次 dispatch，整个弹射循环在单个 shader 内完成。
- **Wavefront**（`--gpu=wavefront`）：同一算法拆成每次弹射 6 个 kernel（raygen / extend / shade / shadow + 间接 dispatch 维护），由常驻的逐像素状态与压缩工作队列连接。

BVH 遍历是在 shader 内实现的，而不是走 `VK_KHR_ray_query`，因为 macOS 上的 MoltenVK 不提供硬件光追扩展。

### 限制

- 只有 `bvh` 这种加速结构会导出 GPU 需要的扁平 BVH；其他 `accelerate` 类型会抛出 `std::runtime_error`。
- 降噪用的特征缓冲（albedo / normal / depth）是在 CPU 上用每像素一条 primary ray 重算的（`compute_gbuffer`），方差估计则来自 GPU 的半缓冲累积——因此 GPU 渲染喂给降噪器的信息与 CPU 路径一致。
- 交互模式目前只做显示，没有"把当前机位写回 XML"的按钮（相机矩阵可由 `CameraController::to_string()` 取得）。

### 独立的 GPU 调试工具

`M_ENABLE_GPU=ON` 时会随 `tiny-renderer` 一起构建（单配置生成器下位于 `build/src/gpu/`，Visual Studio / Xcode 下位于 `build/src/gpu/Release`）：

| 可执行文件 | 用途 |
| --- | --- |
| `gpu-smoketest` | 最小的 Vulkan compute 冒烟测试（不需要场景） |
| `gpu-raytrace-debug <scene.xml> [output.png]` | BVH 遍历 / 法线可视化 |
| `gpu-path-trace <scene.xml> [output.png] [spp]` | 离线 megakernel 路径追踪 |
| `gpu-wavefront <scene.xml> [output.png] [spp]` | 离线 wavefront 路径追踪 |
| `gpu-interactive-verify <scene.xml> [spp]` | 校验"N 帧交互渲染 == 离线 megakernel 的 N spp"，并检查时域滤波是否收敛 |
| `atomic-test` | 设备端原子操作诊断（不需要场景） |



## 场景文件格式

场景是一个以 `<scene>` 为根的 XML 文档。采样器、加速结构、积分器、相机以及可选的降噪器各声明一次，其后是网格与光源。

```xml
<scene>
    <sampler type="independent">
        <integer name="sample_count" value="1024"/>
    </sampler>

    <accelerate type="bvh">
        <integer name="leaf_max" value="5"/>
        <integer name="max_depth" value="100"/>
    </accelerate>

    <integrator type="path">
        <integer name="max_depth" value="5"/>
        <integer name="rr_depth" value="5"/>
    </integrator>

    <camera type="perspective">
        <rfilter type="tent">
            <float name="radius" value="0.5"/>
        </rfilter>
        <transform name="to_world">
            <matrix value="1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1"/>
        </transform>
        <float name="fov" value="35"/>
        <integer name="width" value="1280"/>
        <integer name="height" value="720"/>
    </camera>

    <mesh type="obj">
        <string name="filename" value="models/Mesh001.obj"/>
        <bsdf type="diffuse">
            <texture type="constant">
                <color name="color" value="0.9, 0.9, 0.9"/>
            </texture>
        </bsdf>
    </mesh>

    <emitter type="envmap">
        <texture type="bitmap">
            <string name="filename" value="textures/envmap.hdr"/>
        </texture>
    </emitter>
</scene>
```

完整场景见 `assets/`（`teapot.xml` 是最小的一个）。

### 降噪器

`<denoiser>` 是 `<scene>` 的子元素，渲染结束（CPU 或 GPU）后自动生效。它可以包含一个嵌套的降噪器，后者先运行、作为预滤波——最典型的组合就是先做萤火虫剔除，再做空域滤波：

```xml
<denoiser type="atrous">
    <integer name="iterations" value="5"/>
    <float name="sigma_c" value="4.0"/>
    <denoiser type="outlier"/>
</denoiser>
```

参数如下（括号内为默认值）：

| 类型 | 参数 |
| --- | --- |
| `atrous` | `iterations` (5)、`sigma_c` (4.0，颜色)、`sigma_n` (128.0，法线)、`sigma_d` (1.0，深度)、`sigma_a` (0.05，反照率)、`demodulate` (true)、`threads` (0) |
| `nlm` | `radius` (10)、`patch_radius` (3)、`k` (0.45)、`k_smooth` (1.0)、`alpha` (1.0)、`cross_validate` (true)、`demodulate` (true)、`sigma_n` (0.8)、`sigma_d` (0.6)、`threads` (0) |
| `outlier` | `radius` (1)、`threshold` (2.0)、`threads` (0) |

`threads = 0` 表示使用 `std::thread::hardware_concurrency()`。



## 代码结构

### 基本数据类型

基本数据类型位于`include/core`，主要是模板类，用于定义基本的数据结构以及实现基础的数学库，包括：

- pcg32：随机数生成。
- array：数组，包括向量、点、法向，这三者在进行线性变换时需要区别。
- matrix：矩阵类，包含矩阵的乘法、求逆、转置等常见操作。
- bounding box：包围盒，包含检测、扩展、合并、光线求交等常见操作。
- distribution：对一维概率密度以及二维概率密度进行建模，实现采样、计算概率、重参数化等操作。
- frame：用于世界坐标和局部坐标的变换。
- fresnel：包含导体、介电体的菲涅尔项计算。
- microfacet：微表面模型建模，包含法线分布项、几何遮挡项等。
- intersection：光线与表面交点的数据结构。
- quad：勒让德积分的实现，用于塑料材质次表面散射的建模。
- ray：光线类。
- record：采样数据存储类。
- spectrum：颜色类，目前只实现了 RGB，没有实现光谱。
- tensor：大型张量，存储动态分配，用于材质的存储。
- timer：计时类。
- transform：变换类，包含对向量、点、法向的变换以及透视投影、旋转平移的矩阵构建。
- warp：一系列随机数分布变换函数。
- gpu_scene：`Scene` 展平后可供 GPU 上传的形式——顶点、三角形、BVH 节点、材质、贴图与光源的纯数组，由 `Scene::build_gpu_scene()` 生成。

### 场景组成

可以由 XML 文件指定实例化的类被称之为组件，代码位于`include/components`，用于描述场景组件以及包装路径追踪需要用到的工具。所有类均继承基类`Object`，并通过`ObjectFactory`类注册构造函数，以便于读取 XML 文件时实例化组件类并组成包含关系。下面是组件类的包含关系图：

![Scene](./assets/document/Scene.png)

除了环境光照`Environment map`外，其余光照`Emitter`均依附于几何体`Mesh`，但是场景`Scene`会保留所有光照的共享指针以便于管理光源并做多光源的重要性采样。

### 渲染工具（`include/render`）

- **block**：`ImageBlock`（线程局部累积缓冲，可选携带 AOV 半缓冲）与 `BlockGenerator`（把图像块分发给各线程的任务调度器）。
- **aov**：降噪器需要时，与辐射值一起写入的辅助采样记录（`albedo` / `normal` / `depth`）。
- **framebuffer**：`FrameBufferSet`——从渲染结果流向降噪器的那一组缓冲（颜色、两个半缓冲、方差、特征缓冲）。
- **gbuffer**：在 CPU 上用每像素一条 primary ray 重算几何特征缓冲（albedo、normal、depth），使 GPU 渲染也能获得与 CPU 一致的特征信息。
- **progress**：`ProgressReporter`——`--progress` 的控制台进度条与可选预览窗口。
- **denoise_utils**：降噪器共用的辅助函数（缓冲转换、并行分带调度等）。

### 降噪器（`include/components/denoiser.h`、`src/denoisers`）

`Denoiser` 是 `Object` 的子类，因此和其他组件一样通过 `ObjectFactory` 直接从 XML 实例化。它接收线性辐射空间的渲染结果，返回同样位于线性空间的滤波图像，并可以持有一个嵌套的降噪器作为预滤波。实现有 `atrous`、`nlm`、`outlier` 三种。

### GPU 后端（`include/gpu`、`src/gpu`）

- **vk_context**：Vulkan 实例 / 设备 / 队列的初始化与 shader 模块加载。
- **gpu_buffer**、**gpu_image**：缓冲与 storage image 的封装（降噪器的屏幕大小缓冲用的是 image 而非 buffer）。
- **gpu_scene_upload**：把 `GPUScene` 打包成 shader 绑定的设备缓冲。
- **gpu_renderer**：`--gpu` 使用的离线 GPU 入口（`render_gpu`、`render_gpu_with_aov`）。
- **gpu_session**：`InteractiveSession`——`--interactive` 背后的常驻会话；把所有一次性初始化提到构造函数里，使 `render_frame()` 足够便宜、可以每秒调用 30 次以上。
- **gpu_timings**：交互面板中显示的各 pass GPU 时间戳查询。
- **src/gpu/shaders**：Slang 源码（`path_trace`、`wavefront_*`、`rt_path_trace`、`svgf_*`），构建时交叉编译为 SPIR-V 输出到构建目录。

### GUI（`include/gui`、`src/gui`）

- **preview_window**：`--progress` 使用的被动窗口——只负责显示传进来的图像，不拥有主循环。
- **interactive_window**：**拥有**交互主循环的窗口：采集输入、承载 ImGui 参数面板、上传 GPU 的 RGBA8 输出，并报告某次改动是否使时域累积失效。
- **camera_controller**：与输入无关的相机运动学（orbit / pan / dolly / fly），直接输出可上传的 `camera_to_world` 与 `sample_to_camera` 矩阵。它刻意不去触碰 `Camera` 类，以免把 GLFW 状态带进离线渲染器。



## 多重重要性采样算法

### 介绍

路径追踪主要依赖递归采样来计算渲染方程。渲染方程可以表示为：
$$
L_o(p,\omega_o)=L_e(p,\omega_o)+\int_\Omega f_r(p,\omega_i,\omega_o)L_i(p,\omega_i)|\cos\theta_i|\mathrm{d}\omega_i
$$
路径追踪使用蒙特卡洛积分来采样光路并计算积分。即对于积分$\displaystyle\int f(x)\mathrm{d}\mu(x)$，在积分空间下按照概率分布$X\sim p(x)$采样样本$x_1\cdots x_n$，并通过公式$\displaystyle{\frac{1}{n}\sum_{i=1}^n\frac{f(x_i)}{p(x_i)}}$来近似积分结果。这样计算积分的方差为：$\displaystyle{\frac{1}{n}\left[\int\frac{f(x)^2}{p(x)}\mathrm{d}\mu(x)-(\displaystyle\int f(x)\mathrm{d}\mu(x))^2\right]}$，因此$p(x)$越接近归一化的$f(x)$，$p(x)$方差越小。

然而我们很难得到$f(x)$的解析解，此时有两种采样策略：一种是根据 BSDF 采样，一种是根据光源分布采样，一种常见的做法是将漫反射和光源采样分开：
$$
\begin{align}
L_o(p,\omega_o)
& =L_e(p,\omega_o)+\int_\Omega f_r(p,\omega_i,\omega_o)L_i(p,\omega_i)|\cos\theta_i|\mathrm{d}\omega_i\\
& =L_e(p,\omega_o)+\int_\Omega f_r(p,\omega_i,\omega_o)(L_e(p',\omega_i)+L_o(p'',\omega_i))|\cos\theta_i|\mathrm{d}\omega_i\\
& =L_e(p,\omega_o)+\int_\Omega f_r(p,\omega_i,\omega_o)L_e(p',\omega_i)|\cos\theta_i|\mathrm{d}\omega_i+\\&\int_\Omega f_r(p,\omega_i,\omega_o)L_o(p'',\omega_i)|\cos\theta_i|\mathrm{d}\omega_i\\
\end{align}
$$
然而这样做会导致镜面材质无法采样到光源，究其原因在于镜面反射分布是一个$\delta$分布，连带着积分是一个$\delta$分布的积分，而光源存在一定的面积，很难正好采样在$\delta$分布的点上，从而造成采样光源得到的 BSDF 概率密度为$0$的情况。

多重重要性采样（Multiple Importance Sampling，MIS）的思路是结合多种不同的采样方法，分别采样被积函数的不同部分，并将这些采样点结合起来，以达到接近于最优采样的结果。在MIS 中，为了拟合积分$\displaystyle\int f(x)\mathrm{d}\mu(x)$的结果，我们采用$m$种采样策略，每种采样策略采样$n_i$次，MIS 公式可以表示为：
$$
\begin{align}
& \sum_{i=1}^m\frac{1}{n_i}\sum_{j=1}^{n_i}w_i(x_{i,j})\frac{f(x_{i,j})}{p_i(x_{i,j})}\\
& \text{where }w_i(x_{i,j})=\frac{(n_ip_i(x_{i,j}))^\beta}{\sum_{k=1}^m (n_kp_k(x_{i,j}))^\beta}
\end{align}
$$
那么对于积分$\displaystyle{\int f_r(p,\omega_i,\omega_o)L_e(p',\omega_i)|\cos\theta_i|\mathrm{d}\omega_i}$，可以以$p$的概率用BSDF采样，用$1-p$的概率按照光源采样，这样可以很大程度上缓解这个问题。多种重要性采样可以表示为：
$$
\begin{align}
& L_o(p,\omega_o)=L_e(p,\omega_o)+\hat{L_s}(p,\omega_o)\\
& \hat{L_s}(p,\omega_o)=\left\{
\begin{array}{cl}
& \hat{L_{s1}}(p,\omega_o)=\frac{f_r(p,\omega_i,\omega_o)L_i(p,\omega_i)|\cos\theta_i|}{p(\omega_i)}\text{ if $f_r$ is specular at p}\\
& \hat{L_{s2}}(p,\omega_o)=\hat{E}(p,\omega_o)+\hat{S}(p,\omega_o)\text{ otherwise}
\end{array}
\right.\\
& \hat{E}(p,\omega_o)=w_1\,\hat{E_1}(p,\omega_o)+w_2\,\hat{E_2}(p,\omega_o)\\
& w_1=\frac{p_{emitter}^2(\omega_i)}{p_{emitter}^2(\omega_i)+p_{bsdf}^2(\omega_i)},\quad w_2=\frac{p_{bsdf}^2(\omega_i)}{p_{emitter}^2(\omega_i)+p_{bsdf}^2(\omega_i)}\\
& \hat{E_1}(p,\omega_o)=\frac{f_r(p,\omega_i,\omega_o)L_e(p',\omega_i)|\cos\theta_i|}{p_{emitter}(\omega_i)}\\
& \hat{E_2}(p,\omega_o)=\frac{f_r(p,\omega_i,\omega_o)L_e(p',\omega_i)|\cos\theta_i|}{p_{bsdf}(\omega_i)}\\
& \hat{S}(p,\omega_o)=\frac{f_r(p,\omega_i,\omega_o)L_o(p'',\omega_i)|\cos\theta_i|}{p_{bsdf}(\omega_i)}\\
\end{align}
$$




### 主循环实现

实现位于`src/integrator/path.cpp`的`li`函数，流程如下：

1. 光线初始化：从相机发射一条光线，沿着视图方向出发，与场景进行第一次相交计算，记录交点信息。

   ```cpp
   SurfaceIntersection3f its;
   bool is_intersect = scene->get_accel()->ray_intersect(ray, its, false);
   ```

   its 记录了交点的各种信息，包括：

   ```cpp
   PointType p;  // 交点位置
   Scalar t;  // 光线传播距离
   NormalType n;  // 法向，如果有法向那就根据法向和重心坐标插值，没有就三角形两条边叉乘
   PointType2 uv;  // uv坐标，有uv就根据uv插值，没有那就是三角形的重心坐标
   FrameType shading_frame;  // 法向插值得到的法向
   FrameType geometric_frame;  // 三角形法向
   VectorType wi;  // 局部坐标系下的入射光线
   VectorType dp_du;  // p随u变化率，用于凹凸贴图
   VectorType dp_dv;  // p随u变化率，用于凹凸贴图
   uint32_t primitive_index;  // 三角形id
   std::shared_ptr<Mesh> mesh;  // 指向对应mesh的共享指针
   ```

2. 直接辐射贡献：

   ```cpp
   // ---------------------- Direct emission ----------------------
   
   // If intersect an emitter
   if (is_intersect && (!its.mesh || its.mesh->is_emitter())) {
       DirectionSample3f ds(its, prev_si);
       if (!its.mesh) {
           ds.emitter = scene->get_environment();
       }
       float em_pdf = 0.0f;
   
       if (!prev_bsdf_delta) {
           em_pdf = scene->pdf_emitter_direction(prev_si, ds, valid);
       }
   
       float mis_bsdf = mis_weight(prev_bsdf_pdf, em_pdf);
   
       result += throughput * ds.emitter->eval(its, valid) * mis_bsdf;
   }
   ```

   这里对应公式中的$\hat{E_2}$以及$\hat{L_{s1}}(p,\omega_o)$中的光源部分。

3. 判断是否要继续弹射：

   ```cpp
   bool active_next = depth + 1 < m_max_depth && is_intersect && its.mesh;
   
   if (!active_next) {
   	break;
   }
   ```

   如果达到最大深度或者没有相交点时退出路径追踪过程。

4. 光源采样：

   ```CPP
   std::shared_ptr<BSDF> bsdf = its.mesh->get_bsdf();
   
   // ---------------------- Emitter sampling ----------------------
   bool active_em = bsdf->has_flag(ESmooth);
   
   DirectionSample3f ds;
   Color3f em_weight;
   Vector3f wo;
   
   if (active_em) {
       std::tie(ds, em_weight) = scene->sample_emitter_direction(its, sampler->next2d(), true, active_em);
       active_em &= ds.pdf != 0.0f;
       wo = its.to_local(ds.d);
   }
   ```

   在交点处，根据场景所有发光体的分布，选取一个方向并得到采样权重 。然后将世界坐标的方向向量转换到交点局部坐标系下，以配合接下来的 BSDF 评估。

5. BSDF 评估与 BSDF 采样

   ```cpp
   // ------ Evaluate BSDF * cos(theta) and sample direction -------
   float sample1   = sampler->next1d();
   Point2f sample2 = sampler->next2d();
   
   auto bsdf_val                   = bsdf->eval(its, wo, active_next);
   auto bsdf_pdf                   = bsdf->pdf(its, wo, active_next);
   auto [bsdf_sample, bsdf_weight] = bsdf->sample(its, sample1, sample2, active_next);
   ```

   前两者评估直接采样光源的 BSDF 值以及对应的概率密度，后者为根据 BSDF 分布重要性采样散射光线的方向、BSDF 值以及概率密度。

6. 多重重要性采样 (MIS) 整合

   ```cpp
   // --------------- Emitter sampling contribution ----------------
   if (active_em) {
       float mis_em = ds.delta ? 1.0f : mis_weight(ds.pdf, bsdf_pdf);
       result += throughput * bsdf_val * em_weight * mis_em;
   }
   ```

   这里对应公式中的$\hat{E_1}$中的光源部分。

7. BSDF 采样得到新方向并更新吞吐量

   ```cpp
   // ---------------------- BSDF sampling ----------------------
   ray = its.spawn_ray(its.to_world(bsdf_sample.wo));
   
   // ------ Update loop variables based on current interaction ------
   throughput *= bsdf_weight;
   eta *= bsdf_sample.eta;
   valid_ray |= valid && its.is_valid();
   
   // Information about the current vertex needed by the next iteration
   prev_si         = its;
   prev_bsdf_pdf   = bsdf_sample.pdf;
   prev_bsdf_delta = bsdf_sample.delta;
   ```

   这里对应公式中的$\hat{S}$以及$\hat{L_{s1}}(p,\omega_o)$中的散射部分，通过循环采样的方式累加。

8. 俄罗斯轮盘赌 (Russian Roulette)

   ```cpp
   // -------------------- Stopping criterion ---------------------
   depth += 1;
   float throughput_max = throughput.max_value();
   float rr_prob        = M_MIN(throughput_max * eta * eta, 0.95f);
   bool rr_active       = depth >= m_rr_depth;
   bool rr_continue     = sampler->next1d() < rr_prob;
   
   if (rr_active) {
   	throughput *= 1.0f / rr_prob;
   }
   
   valid = (!rr_active || rr_continue) && throughput_max != 0.0f;
   ```

   它的作用是在路径足够长时，随机截断一些能量很弱的路径，从而避免无限递归和无意义计算，同时保证结果无偏。具体而言，其在 `depth >= m_rr_depth`（比如第 5、6 次弹射以后）启用俄罗斯轮盘赌；计算一个继续概率 `rr_prob = min(throughput_max * eta^2, 0.95)`，如果随机数 `sampler->next1d()` 小于这个 `rr_prob`，才继续弹射，否则就终止。若继续了，则要把 `throughput` 再乘上 `1 / rr_prob`，补偿这一随机过程带来的期望下降，保证最终结果不带偏差。



### 光源的主动采样策略

对于场景多光源的情况，首先以均匀分布随机采样一个多面体光源。在多面体光源上，将三角形面积与三角形三顶点平均luminance之积作为概率密度函数进行重要性采样，从而对场景中的光源进行重要性采样。具体实现见`scene.cpp`的`sample_emitter_direction`函数与`area.cpp`的实现。



## 结果展示

<img src="./assets/box/box.png" alt="box"/>

<img src="./assets/ball/ball.png" alt="ball"/>

<img src="./assets/mis/mis.png" alt="mis"/>

<img src="./assets/bathroom2/bathroom2.png" alt="bathroom2"/>

<img src="./assets/living-room/living-room.png" alt="living-room"/>

<img src="./assets/teapot/teapot.png" alt="teapot"/>

<img src="./assets/bidir/bidir.png" alt="bidir"/>

### 降噪

每个场景目录下还保存了未降噪的渲染结果（`*_noisy.png`，只要运行了降噪器就会生成），可以直接 A/B 对比降噪效果：

| 未降噪（`*_noisy.png`） | 降噪后（`*.png`） |
| --- | --- |
| <img src="./assets/living-room/living-room_noisy.png" alt="living-room 未降噪"/> | <img src="./assets/living-room/living-room.png" alt="living-room 降噪后"/> |
| <img src="./assets/bathroom2/bathroom2_noisy.png" alt="bathroom2 未降噪"/> | <img src="./assets/bathroom2/bathroom2.png" alt="bathroom2 降噪后"/> |
| <img src="./assets/teapot/teapot_noisy.png" alt="teapot 未降噪"/> | <img src="./assets/teapot/teapot.png" alt="teapot 降噪后"/> |



## 未来更新

1. 体渲染：散射介质、相函数以及对应的路径追踪器
2. 更多积分器：双向路径追踪、Metropolis Light Transport
3. 更多材质：法向贴图、BSSRDF（[Position-Free Monte Carlo Simulation for Arbitrary Layered BSDFs](https://shuangz.com/projects/layered-sa18/)）、毛发等
4. 更多类型的几何体，曲面曲线的实现
5. ~~可交互GUI~~ 已完成——`--interactive`，见[实时交互模式](#实时交互模式)。
6. ~~光线追踪降噪~~ 已完成——离线的 `atrous` / `nlm` / `outlier` 降噪器，以及交互模式中的 SVGF。
7. GPU 后端：在支持的平台上改用硬件光追（`VK_KHR_ray_query`），以及与 CPU 路径完全对齐的材质/体积支持。
8. 把交互模式中找到的相机机位写回场景 XML。



## 未修复Bug

1. 图像分辨率不是$2^x$时保存图像会出错。（历史问题：现在 `Bitmap::save_png()` 通过 `stbi_write_png` 写出，没有 2 的幂限制——该问题可能已不复现，需要重新验证。）
2. 双层材质`Smooth`标签的识别。
3. 缺少 GPU 后端或预览窗口时 `--interactive` 不可用；若 `glfw3`/OpenGL 缺失，CMake 只会打印一条状态信息并继续，因此问题要到运行时才暴露（Windows/Apple 会退回编译自带的 `ext/glfw`，不会走到这一步）。
4. GPU 后端只支持能导出扁平 BVH 的加速结构（目前只有 `bvh`），其他类型会抛出 `std::runtime_error`。
5. 命令行 usage 中仍然列着 `--no-gui`，但它并未实现；预览窗口是由 `--progress` 打开的。
6. 交互模式不输出图像，也还没有"把当前机位复制回 XML"的按钮。
