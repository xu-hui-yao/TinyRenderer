#include <components/bsdf.h>
#include <components/object.h>

M_NAMESPACE_BEGIN

Object::EClassType BSDF::get_class_type() const { return EBSDF; }

bool BSDF::has_flag(BSDFFlags bsdf_flags) const { return (m_flags & bsdf_flags) != 0; }

BSDFFlags BSDF::get_flag() const { return m_flags; }

M_NAMESPACE_END