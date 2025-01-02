#include <components/sampler.h>
#include <core/pcg32.h>

M_NAMESPACE_BEGIN
/**
 * Independent sampling - returns independent uniformly distributed
 * random numbers on <tt>[0, 1)x[0, 1)</tt>.
 *
 * This class is essentially just a wrapper around the pcg32 pseudorandom
 * number generator. For more details on what sample generators do in
 * general, refer to the \ref Sampler class.
 */
class Independent : public Sampler {
public:
    explicit Independent(const PropertyList &property_list) {
        m_sample_count = static_cast<size_t>(property_list.get_integer("sample_count", 1));
        this->m_name   = property_list.get_string("name", "independent");
    }

    void construct() override {
#ifdef M_DEBUG
        std::cout << "Construct " << class_type_name(get_class_type()) << std::endl;
#endif
    }

    void add_child(const std::shared_ptr<Object> &child) override {}

    ~Independent() override = default;

    [[nodiscard]] std::shared_ptr<Sampler> clone() const override {
        std::shared_ptr<Independent> cloned(new Independent());
        cloned->m_sample_count = m_sample_count;
        cloned->m_random       = m_random;
        return cloned;
    }

    void prepare(const ImageBlock &block) override { m_random.seed(block.get_offset().x(), block.get_offset().y()); }

    void generate() override { /* No-op for this sampler */ }

    void advance() override { /* No-op for this sampler */ }

    float next1d() override { return m_random.next_float(); }

    Point2f next2d() override { return { m_random.next_float(), m_random.next_float() }; }

    [[nodiscard]] std::string to_string() const override {
        return std::string("Independent[sample_count=") + std::to_string(m_sample_count) + std::string("]");
    }

protected:
    Independent() {}

private:
    pcg32 m_random;
};

REGISTER_CLASS(Independent, "independent");

M_NAMESPACE_END
