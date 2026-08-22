#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert Mitsuba 3 XML scene files (assets/mitsuba3/<scene>/scene.xml) into
the TinyRenderer XML dialect (assets/mitsuba3/<scene>/<scene>.xml).

This is a one-off migration helper script; not part of the shipped codebase.
"""
import copy
import os
import re
import sys
import xml.etree.ElementTree as ET

NOTES = []  # collected (scene_name, message) incompatibility notes


def note(scene, msg):
    NOTES.append((scene, msg))
    print(f"[NOTE][{scene}] {msg}")


# ---------------------------------------------------------------------------
# Low level helpers
# ---------------------------------------------------------------------------

def mk(tag, **attrib):
    e = ET.Element(tag)
    for k, v in attrib.items():
        e.set(k, str(v))
    return e


def mk_const_color_texture(value_str):
    t = mk('texture', type='constant')
    t.append(mk('color', name='color', value=value_str))
    return t


def xyz_to_value(node, defx='0', defy='0', defz='0'):
    x = node.get('x', defx)
    y = node.get('y', defy)
    z = node.get('z', defz)
    return f"{x}, {y}, {z}"


def copy_property(node, scene):
    """Generic recursive copy of a *property* (non object) XML node,
    translating mitsuba-specific spellings into tinyrenderer's dialect."""
    tag = node.tag

    if tag == 'rgb':
        new = mk('color')
        if node.get('name') is not None:
            new.set('name', node.get('name'))
        new.set('value', node.get('value'))
        return new

    if tag == 'spectrum':
        # Not used by any of our scenes but handle gracefully: treat the
        # (uniform) value as a grayscale color if present.
        val = node.get('value', '1')
        new = mk('color')
        if node.get('name') is not None:
            new.set('name', node.get('name'))
        new.set('value', f"{val}, {val}, {val}")
        return new

    if tag in ('point', 'vector'):
        new = mk(tag)
        if node.get('name') is not None:
            new.set('name', node.get('name'))
        if 'value' in node.attrib:
            new.set('value', node.get('value'))
        else:
            new.set('value', xyz_to_value(node))
        return new

    if tag == 'rotate':
        new = mk('rotate')
        new.set('angle', node.get('angle'))
        new.set('axis', xyz_to_value(node))
        return new

    if tag == 'scale':
        new = mk('scale')
        if 'value' in node.attrib:
            new.set('value', node.get('value'))
        else:
            new.set('value', xyz_to_value(node, '1', '1', '1'))
        return new

    if tag == 'translate':
        new = mk('translate')
        if 'value' in node.attrib:
            new.set('value', node.get('value'))
        else:
            new.set('value', xyz_to_value(node))
        return new

    if tag == 'lookat' or tag == 'look_at':
        new = mk('look_at')
        for a in ('origin', 'target', 'up'):
            if a in node.attrib:
                new.set(a, node.get(a))
        return new

    if tag == 'transform':
        new = mk('transform')
        if node.get('name') is not None:
            new.set('name', node.get('name'))
        for c in node:
            new.append(copy_property(c, scene))
        return new

    if tag == 'texture':
        return convert_texture(node, scene)

    # float / integer / boolean / string / matrix -> passthrough clone
    new = ET.Element(tag, attrib=dict(node.attrib))
    for c in node:
        new.append(copy_property(c, scene))
    return new


def convert_texture(node, scene):
    ttype = node.get('type')
    if ttype == 'bitmap':
        new = mk('texture', type='bitmap')
        for c in node:
            new.append(copy_property(c, scene))
        return new

    if ttype == 'checkerboard':
        new = mk('texture', type='checkerboard')
        for c in node:
            if c.tag == 'transform' and c.get('name') == 'to_uv':
                scale_el = None
                for cc in c:
                    if cc.tag == 'scale':
                        scale_el = cc
                        break
                if scale_el is not None:
                    su = scale_el.get('x', scale_el.get('value', '1'))
                    sv = scale_el.get('y', scale_el.get('value', '1'))
                    new.append(mk('float', name='scale_u', value=su))
                    new.append(mk('float', name='scale_v', value=sv))
                note(scene, "checkerboard <transform name=\"to_uv\"> converted to scale_u/scale_v "
                             "(rotation/shear part of the transform, if any, is NOT supported)")
            else:
                new.append(copy_property(c, scene))
        return new

    if ttype == 'constant':
        new = mk('texture', type='constant')
        for c in node:
            new.append(copy_property(c, scene))
        return new

    if ttype == 'volume':
        note(scene, "texture type=\"volume\" (gridvolume-backed texture) is NOT supported; "
                     "substituted with a flat gray constant texture")
        return mk_const_color_texture("0.5, 0.5, 0.5")

    note(scene, f"texture type=\"{ttype}\" is NOT supported; substituted with a flat gray constant texture")
    return mk_const_color_texture("0.5, 0.5, 0.5")


def prop_to_texture(node, scene):
    """Turn a mitsuba property node (float / rgb / texture) that represents
    a color/scalar into a tinyrenderer <texture> child."""
    if node.tag == 'texture':
        t = convert_texture(node, scene)
        if 'name' in t.attrib:
            del t.attrib['name']
        return t
    if node.tag == 'rgb':
        return mk_const_color_texture(node.get('value'))
    if node.tag == 'spectrum':
        v = node.get('value', '1')
        return mk_const_color_texture(f"{v}, {v}, {v}")
    if node.tag == 'float':
        v = node.get('value')
        return mk_const_color_texture(f"{v}, {v}, {v}")
    raise ValueError(f"Cannot convert <{node.tag}> to a texture")


def find_prop(children, name):
    for c in children:
        if c.get('name') == name:
            return c
    return None


def find_bsdf_children(children):
    return [c for c in children if c.tag == 'bsdf']


# ---------------------------------------------------------------------------
# BSDF conversion
# ---------------------------------------------------------------------------

def convert_bsdf(elem, scene):
    btype = elem.get('type')
    children = list(elem)
    new = ET.Element('bsdf')

    if btype == 'twosided':
        new.set('type', 'two_sided')
        for b in find_bsdf_children(children)[:2]:
            new.append(convert_bsdf(b, scene))
        return new

    if btype == 'mask':
        new.set('type', 'mask')
        opacity = find_prop(children, 'opacity')
        if opacity is None:
            note(scene, "mask bsdf missing 'opacity' property; defaulted to fully opaque (1.0)")
            new.append(mk_const_color_texture("1, 1, 1"))
        else:
            new.append(prop_to_texture(opacity, scene))
        nested = find_bsdf_children(children)
        if nested:
            new.append(convert_bsdf(nested[0], scene))
        return new

    if btype == 'bumpmap':
        new.set('type', 'bumpmap')
        map_prop = find_prop(children, 'map')
        if map_prop is not None:
            tex = convert_texture(map_prop, scene) if map_prop.tag == 'texture' else prop_to_texture(map_prop, scene)
            if 'name' in tex.attrib:
                del tex.attrib['name']
            new.append(tex)
        nested = find_bsdf_children(children)
        if nested:
            new.append(convert_bsdf(nested[0], scene))
        scale_prop = find_prop(children, 'scale')
        if scale_prop is not None:
            new.append(copy_property(scale_prop, scene))
        return new

    if btype == 'diffuse':
        new.set('type', 'diffuse')
        refl = find_prop(children, 'reflectance')
        if refl is None:
            refl = next((c for c in children if c.tag in ('texture', 'rgb', 'float', 'spectrum')), None)
        if refl is None:
            note(scene, "diffuse bsdf missing reflectance; defaulted to 0.5 gray")
            new.append(mk_const_color_texture("0.5, 0.5, 0.5"))
        else:
            new.append(prop_to_texture(refl, scene))
        return new

    if btype == 'conductor':
        new.set('type', 'conductor')
        material = find_prop(children, 'material')
        eta = find_prop(children, 'eta')
        k = find_prop(children, 'k')
        specr = find_prop(children, 'specular_reflectance')
        if material is not None:
            mat_val = material.get('value', 'none')
            if mat_val != 'none':
                note(scene, f"conductor material=\"{mat_val}\" preset is NOT supported "
                             "(only explicit eta/k are supported); substituted with a perfect mirror (eta=0, k=1)")
            eta_tex = mk_const_color_texture("0, 0, 0")
            k_tex = mk_const_color_texture("1, 1, 1")
        else:
            eta_tex = prop_to_texture(eta, scene) if eta is not None else mk_const_color_texture("0, 0, 0")
            k_tex = prop_to_texture(k, scene) if k is not None else mk_const_color_texture("1, 1, 1")
        new.append(eta_tex)
        new.append(k_tex)
        if specr is not None:
            new.append(prop_to_texture(specr, scene))
        return new

    if btype == 'roughconductor':
        new.set('type', 'roughconductor')
        material = find_prop(children, 'material')
        eta = find_prop(children, 'eta')
        k = find_prop(children, 'k')
        alpha = find_prop(children, 'alpha')
        specr = find_prop(children, 'specular_reflectance')
        distribution = find_prop(children, 'distribution')
        if distribution is not None and distribution.get('value') != 'ggx':
            note(scene, f"roughconductor distribution=\"{distribution.get('value')}\" is NOT supported "
                         "(only a GGX-like distribution is implemented); ignored")
        if material is not None:
            mat_val = material.get('value', 'none')
            if mat_val != 'none':
                note(scene, f"roughconductor material=\"{mat_val}\" preset is NOT supported; "
                             "substituted with a perfect mirror (eta=0, k=1)")
            eta_tex = mk_const_color_texture("0, 0, 0")
            k_tex = mk_const_color_texture("1, 1, 1")
        else:
            eta_tex = prop_to_texture(eta, scene) if eta is not None else mk_const_color_texture("0, 0, 0")
            k_tex = prop_to_texture(k, scene) if k is not None else mk_const_color_texture("1, 1, 1")
        alpha_tex = prop_to_texture(alpha, scene) if alpha is not None else mk_const_color_texture("0.1, 0.1, 0.1")
        new.append(eta_tex)
        new.append(k_tex)
        new.append(alpha_tex)
        if specr is not None:
            new.append(prop_to_texture(specr, scene))
        return new

    if btype in ('dielectric', 'thindielectric'):
        new.set('type', btype)
        for key in ('int_ior', 'ext_ior'):
            p = find_prop(children, key)
            if p is not None:
                new.append(copy_property(p, scene))
        specr = find_prop(children, 'specular_reflectance')
        spect = find_prop(children, 'specular_transmittance')
        if specr is not None:
            new.append(prop_to_texture(specr, scene))
        if spect is not None:
            new.append(prop_to_texture(spect, scene))
        return new

    if btype == 'roughdielectric':
        new.set('type', 'roughdielectric')
        for key in ('int_ior', 'ext_ior'):
            p = find_prop(children, key)
            if p is not None:
                new.append(copy_property(p, scene))
        distribution = find_prop(children, 'distribution')
        if distribution is not None and distribution.get('value') != 'ggx':
            note(scene, f"roughdielectric distribution=\"{distribution.get('value')}\" is NOT supported; ignored")
        alpha = find_prop(children, 'alpha')
        new.append(prop_to_texture(alpha, scene) if alpha is not None else mk_const_color_texture("0.1, 0.1, 0.1"))
        specr = find_prop(children, 'specular_reflectance')
        spect = find_prop(children, 'specular_transmittance')
        if specr is not None:
            new.append(prop_to_texture(specr, scene))
        if spect is not None:
            new.append(prop_to_texture(spect, scene))
        return new

    if btype == 'plastic':
        new.set('type', 'plastic')
        for key in ('int_ior', 'ext_ior', 'nonlinear'):
            p = find_prop(children, key)
            if p is not None:
                new.append(copy_property(p, scene))
        diffr = find_prop(children, 'diffuse_reflectance')
        new.append(prop_to_texture(diffr, scene) if diffr is not None else mk_const_color_texture("0.5, 0.5, 0.5"))
        specr = find_prop(children, 'specular_reflectance')
        if specr is not None:
            new.append(prop_to_texture(specr, scene))
        return new

    if btype == 'roughplastic':
        new.set('type', 'roughplastic')
        distribution = find_prop(children, 'distribution')
        if distribution is not None and distribution.get('value') != 'ggx':
            note(scene, f"roughplastic distribution=\"{distribution.get('value')}\" is NOT supported; ignored")
        for key in ('int_ior', 'ext_ior', 'nonlinear', 'alpha'):
            p = find_prop(children, key)
            if p is not None:
                new.append(copy_property(p, scene))
        diffr = find_prop(children, 'diffuse_reflectance')
        new.append(prop_to_texture(diffr, scene) if diffr is not None else mk_const_color_texture("0.5, 0.5, 0.5"))
        specr = find_prop(children, 'specular_reflectance')
        if specr is not None:
            new.append(prop_to_texture(specr, scene))
        return new

    note(scene, f"bsdf type=\"{btype}\" is NOT supported; substituted with a plain diffuse gray material")
    fallback = mk('bsdf', type='diffuse')
    fallback.append(mk_const_color_texture("0.5, 0.5, 0.5"))
    return fallback


# ---------------------------------------------------------------------------
# Emitter conversion
# ---------------------------------------------------------------------------

def convert_emitter(elem, scene):
    etype = elem.get('type')
    children = list(elem)
    if etype == 'area':
        new = mk('emitter', type='area')
        radiance = find_prop(children, 'radiance')
        if radiance is None:
            radiance = next((c for c in children if c.tag in ('texture', 'rgb', 'float', 'spectrum')), None)
        new.append(prop_to_texture(radiance, scene) if radiance is not None else mk_const_color_texture("1, 1, 1"))
        return new

    if etype == 'envmap':
        new = mk('emitter', type='envmap')
        filename = find_prop(children, 'filename')
        for c in children:
            if c.tag == 'transform':
                new.append(copy_property(c, scene))
        if filename is not None:
            tex = mk('texture', type='bitmap')
            tex.append(mk('string', name='filename', value=filename.get('value')))
            new.append(tex)
        scale = find_prop(children, 'scale')
        if scale is not None:
            note(scene, "envmap 'scale' property is NOT supported; ignored (bake the scale into the "
                         "texture/exposure manually if needed)")
        return new

    note(scene, f"emitter type=\"{etype}\" is NOT supported and was DROPPED from the scene")
    return None


# ---------------------------------------------------------------------------
# Shape / mesh conversion
# ---------------------------------------------------------------------------

def convert_shape(elem, scene, bsdf_defs, texture_defs, ply_map):
    stype = elem.get('type')
    children = list(elem)

    if stype == 'ply':
        new = mk('mesh', type='obj')
    elif stype == 'disk':
        new = mk('mesh', type='rectangle')
        note(scene, "shape type=\"disk\" is NOT supported; substituted with a flat rectangle "
                     "of the same transform (visual shape differs: square vs. circle)")
    elif stype in ('obj', 'rectangle', 'sphere', 'cube'):
        new = mk('mesh', type=stype)
    else:
        note(scene, f"shape type=\"{stype}\" is NOT supported; shape was DROPPED from the scene")
        return None

    for c in children:
        if c.tag == 'ref':
            target_id = c.get('id')
            if target_id in bsdf_defs:
                new.append(convert_bsdf(bsdf_defs[target_id], scene))
            elif target_id in texture_defs:
                t = convert_texture(texture_defs[target_id], scene)
                if 'name' in t.attrib:
                    del t.attrib['name']
                new.append(t)
            else:
                note(scene, f"<ref id=\"{target_id}\"> could not be resolved; ignored")
        elif c.tag == 'bsdf':
            new.append(convert_bsdf(c, scene))
        elif c.tag == 'emitter':
            em = convert_emitter(c, scene)
            if em is not None:
                new.append(em)
        elif c.tag == 'medium':
            note(scene, f"participating <medium type=\"{c.get('type')}\"> is NOT supported; dropped "
                         "(affected object will render as clear/solid instead of a scattering medium)")
            continue
        elif c.tag == 'string' and c.get('name') == 'filename':
            val = c.get('value')
            if stype == 'ply':
                val = ply_map.get(val, val)
            new.append(mk('string', name='filename', value=val))
        else:
            new.append(copy_property(c, scene))
    return new


# ---------------------------------------------------------------------------
# Top level scene conversion
# ---------------------------------------------------------------------------

def substitute_defaults(text):
    defaults = {}
    for m in re.finditer(r'<default\s+name="([^"]+)"\s+value="([^"]*)"\s*/>', text):
        defaults[m.group(1)] = m.group(2)

    def repl(match):
        key = match.group(1)
        return defaults.get(key, match.group(0))

    text = re.sub(r'\$([A-Za-z_][A-Za-z0-9_]*)', repl, text)
    text = re.sub(r'<default\s+[^/]*/>\s*', '', text)
    return text, defaults


def convert_scene(src_path, dst_path, scene_name):
    with open(src_path, 'r', encoding='utf-8') as f:
        text = f.read()

    text, defaults = substitute_defaults(text)
    root = ET.fromstring(text)

    bsdf_defs = {}
    texture_defs = {}
    shapes = []
    top_emitters = []
    sensor_elem = None
    integrator_elem = None

    def collect_ids(node):
        for c in node.iter():
            if c is node:
                continue
            if c.tag == 'bsdf' and c.get('id') and c.get('id') not in bsdf_defs:
                bsdf_defs[c.get('id')] = c
            elif c.tag == 'texture' and c.get('id') and c.get('id') not in texture_defs:
                texture_defs[c.get('id')] = c

    collect_ids(root)

    for c in list(root):
        if c.tag == 'bsdf':
            # Top-level bsdf definitions are only used via <ref id="..."/>
            # elsewhere; ids (including those nested inside unreferenced
            # wrapper bsdfs, e.g. a top-level bumpmap with no id whose
            # child *does* have one) were already collected above.
            continue
        elif c.tag == 'texture':
            continue
        elif c.tag == 'shape':
            shapes.append(c)
        elif c.tag == 'emitter':
            top_emitters.append(c)
        elif c.tag == 'sensor':
            sensor_elem = c
        elif c.tag == 'integrator':
            integrator_elem = c
        # 'default' already stripped

    ply_map = {}
    if any(s.get('type') == 'ply' for s in shapes):
        ply_map = convert_ply_meshes(os.path.dirname(src_path), scene_name)

    out = mk('scene')

    # sampler + film pulled out of <sensor>
    sampler_elem = None
    film_elem = None
    if sensor_elem is not None:
        for c in sensor_elem:
            if c.tag == 'sampler':
                sampler_elem = c
            elif c.tag == 'film':
                film_elem = c

    if sampler_elem is not None:
        new_sampler = mk('sampler', type=sampler_elem.get('type', 'independent'))
        for c in sampler_elem:
            new_sampler.append(copy_property(c, scene_name))
        out.append(new_sampler)
    else:
        s = mk('sampler', type='independent')
        s.append(mk('integer', name='sample_count', value='64'))
        out.append(s)

    accel = mk('accelerate', type='bvh')
    accel.append(mk('integer', name='leaf_max', value='8'))
    accel.append(mk('integer', name='max_depth', value='64'))
    out.append(accel)

    if integrator_elem is not None:
        itype = integrator_elem.get('type', 'path')
        if itype != 'path':
            note(scene_name, f"integrator type=\"{itype}\" is NOT supported (only 'path' is implemented, "
                              "participating-media transport is not simulated); substituted with 'path'")
            itype = 'path'
        new_integrator = mk('integrator', type=itype)
        for c in integrator_elem:
            if c.get('name') == 'hide_emitters':
                note(scene_name, "integrator property 'hide_emitters' is NOT supported; ignored "
                                  "(camera-visible emitters will be rendered instead of hidden)")
                continue
            new_integrator.append(copy_property(c, scene_name))
        out.append(new_integrator)
    else:
        i = mk('integrator', type='path')
        i.append(mk('integer', name='max_depth', value='65'))
        out.append(i)

    camera = mk('camera', type=sensor_elem.get('type', 'perspective') if sensor_elem is not None else 'perspective')
    if film_elem is not None:
        rfilter_elem = None
        width_elem = None
        height_elem = None
        for c in film_elem:
            if c.tag == 'rfilter':
                rfilter_elem = c
            elif c.tag == 'integer' and c.get('name') == 'width':
                width_elem = c
            elif c.tag == 'integer' and c.get('name') == 'height':
                height_elem = c
        if rfilter_elem is not None:
            rtype = rfilter_elem.get('type', 'tent')
            if rtype != 'tent':
                note(scene_name, f"rfilter type=\"{rtype}\" is NOT supported (only 'tent' is implemented); "
                                  "substituted with 'tent'")
                rtype = 'tent'
            new_rfilter = mk('rfilter', type=rtype)
            radius_prop = find_prop(list(rfilter_elem), 'radius')
            new_rfilter.append(mk('float', name='radius', value=radius_prop.get('value') if radius_prop is not None else '1.0'))
            camera.append(new_rfilter)
        else:
            rf = mk('rfilter', type='tent')
            rf.append(mk('float', name='radius', value='1.0'))
            camera.append(rf)
    else:
        rf = mk('rfilter', type='tent')
        rf.append(mk('float', name='radius', value='1.0'))
        camera.append(rf)

    if sensor_elem is not None:
        for c in sensor_elem:
            if c.tag in ('sampler', 'film'):
                continue
            camera.append(copy_property(c, scene_name))
    if film_elem is not None:
        if width_elem is not None:
            camera.append(mk('integer', name='width', value=width_elem.get('value')))
        if height_elem is not None:
            camera.append(mk('integer', name='height', value=height_elem.get('value')))

    out.append(camera)

    for s in shapes:
        m = convert_shape(s, scene_name, bsdf_defs, texture_defs, ply_map)
        if m is not None:
            out.append(m)

    for e in top_emitters:
        conv = convert_emitter(e, scene_name)
        if conv is not None:
            out.append(conv)

    indent(out)
    tree = ET.ElementTree(out)
    tree.write(dst_path, encoding='utf-8', xml_declaration=False)
    print(f"Wrote {dst_path}")


def indent(elem, level=0):
    i = "\n" + level * "    "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "    "
        for child in elem:
            indent(child, level + 1)
        if not child.tail or not child.tail.strip():
            child.tail = i
    if level and (not elem.tail or not elem.tail.strip()):
        elem.tail = i


# ---------------------------------------------------------------------------
# PLY -> OBJ conversion (binary little-endian PLY, as used by the lego scene)
# ---------------------------------------------------------------------------

def convert_ply_meshes(scene_dir, scene_name):
    import struct
    mesh_dir = os.path.join(scene_dir, 'meshes')
    if not os.path.isdir(mesh_dir):
        return {}
    out_dir = os.path.join(scene_dir, 'meshes_obj')
    os.makedirs(out_dir, exist_ok=True)

    ply_map = {}
    files = sorted(f for f in os.listdir(mesh_dir) if f.lower().endswith('.ply'))
    note(scene_name, f"converted {len(files)} binary PLY meshes to OBJ (shape type=\"ply\" is NOT supported "
                       "by the renderer, which only reads Wavefront OBJ)")
    for fname in files:
        src = os.path.join(mesh_dir, fname)
        dst_name = os.path.splitext(fname)[0] + '.obj'
        dst = os.path.join(out_dir, dst_name)
        with open(src, 'rb') as f:
            data = f.read()

        header_end = data.find(b'end_header\n') + len(b'end_header\n')
        header = data[:header_end].decode('ascii', errors='ignore')
        body = data[header_end:]

        fmt_line = re.search(r'format\s+(\S+)', header)
        fmt = fmt_line.group(1) if fmt_line else 'binary_little_endian'
        endian = '<' if 'little' in fmt else '>'

        n_vertex = int(re.search(r'element vertex (\d+)', header).group(1))
        n_face = int(re.search(r'element face (\d+)', header).group(1))

        vprops = re.findall(r'property\s+(\S+)\s+(\S+)\n(?=(?:property|element|end_header))', header + 'element\n')
        # Fallback simpler parse: capture all "property <type> <name>" lines that occur
        # before the 'element face' section.
        vertex_section = header.split('element vertex', 1)[1].split('element face', 1)[0]
        vprop_list = re.findall(r'property\s+(\S+)\s+(\S+)', vertex_section)

        type_size = {'float': 4, 'float32': 4, 'double': 8, 'uchar': 1, 'uint8': 1,
                     'int': 4, 'int32': 4, 'uint': 4, 'uint32': 4, 'short': 2, 'ushort': 2}
        vertex_stride = sum(type_size.get(t, 4) for t, _ in vprop_list)

        positions = []
        normals = []
        uvs = []
        offset = 0
        for _ in range(n_vertex):
            values = {}
            base = offset
            cur = 0
            for t, pname in vprop_list:
                size = type_size.get(t, 4)
                raw = body[base + cur: base + cur + size]
                if t in ('float', 'float32'):
                    values[pname] = struct.unpack_from(endian + 'f', raw)[0]
                elif t in ('double',):
                    values[pname] = struct.unpack_from(endian + 'd', raw)[0]
                else:
                    values[pname] = int.from_bytes(raw, 'little' if endian == '<' else 'big')
                cur += size
            offset += vertex_stride
            positions.append((values.get('x', 0.0), values.get('y', 0.0), values.get('z', 0.0)))
            if 'nx' in values:
                normals.append((values.get('nx', 0.0), values.get('ny', 0.0), values.get('nz', 0.0)))
            if 'u' in values or 's' in values:
                uvs.append((values.get('u', values.get('s', 0.0)), values.get('v', values.get('t', 0.0))))

        faces = []
        for _ in range(n_face):
            count = body[offset]
            offset += 1
            idx = struct.unpack_from(endian + f'{count}i', body, offset)
            offset += 4 * count
            if count == 3:
                faces.append(idx)
            elif count == 4:
                faces.append((idx[0], idx[1], idx[2]))
                faces.append((idx[0], idx[2], idx[3]))

        with open(dst, 'w') as out:
            for p in positions:
                out.write(f"v {p[0]} {p[1]} {p[2]}\n")
            for uv in uvs:
                out.write(f"vt {uv[0]} {uv[1]}\n")
            for n in normals:
                out.write(f"vn {n[0]} {n[1]} {n[2]}\n")
            has_uv = len(uvs) == len(positions)
            has_n = len(normals) == len(positions)
            for tri in faces:
                parts = []
                for vi in tri:
                    s = str(vi + 1)
                    if has_uv or has_n:
                        s += "/" + (str(vi + 1) if has_uv else "")
                    if has_n:
                        s += "/" + str(vi + 1)
                    parts.append(s)
                out.write("f " + " ".join(parts) + "\n")

        ply_map[f"meshes/{fname}"] = f"meshes_obj/{dst_name}"
    return ply_map


if __name__ == '__main__':
    base = '/Users/xuhuiyao/Desktop/workfiled/TinyRenderer/assets/mitsuba3'
    scenes = sys.argv[1:] if len(sys.argv) > 1 else sorted(
        d for d in os.listdir(base) if os.path.isdir(os.path.join(base, d)))
    for scene in scenes:
        src = os.path.join(base, scene, 'scene.xml')
        if not os.path.isfile(src):
            continue
        dst = os.path.join(base, scene, f'{scene}.xml')
        print(f"=== Converting {scene} ===")
        try:
            convert_scene(src, dst, scene)
        except Exception as e:
            print(f"[ERROR] Failed to convert {scene}: {e}")
            import traceback
            traceback.print_exc()

    print("\n\n===== SUMMARY OF INCOMPATIBILITIES =====")
    for scene, msg in NOTES:
        print(f"- [{scene}] {msg}")
