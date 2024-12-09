#pragma once

#include <components/object.h>

M_NAMESPACE_BEGIN

std::shared_ptr<Object> load_from_xml(const std::string& filename);

M_NAMESPACE_END
