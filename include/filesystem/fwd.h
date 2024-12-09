#pragma once

#if !defined(NAMESPACE_BEGIN)
#define NAMESPACE_BEGIN(name) namespace name {
#endif
#if !defined(NAMESPACE_END)
#define NAMESPACE_END(name) }
#endif

NAMESPACE_BEGIN(filesystem)

class path;
class resolver;

NAMESPACE_END(filesystem)
