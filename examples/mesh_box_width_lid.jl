using CompScienceMeshes


function mesh_box_width_lid(size, h1, h2)
    width, height, depth = size, size, size

    m = -meshrectangle(width, height, h1, 3)
    m_front = rotate(m, 0.5*pi*[1,0,0])
    m_left  = rotate(m, -0.5*pi*[0,1,0])
    m_back  = -translate(m_front, [0,height,0])
    m_right = -translate(m_left, [width,0,0])
    m_bottom = -m
    m_open = weld(m_front, m_left, m_back, m_right, m_bottom)

    m_fine = -meshrectangle(width, height, h2)
    m_top  = translate(m_fine, [0,0,height])

    # m_coarse = -meshrectangle(width, height, h1)
    # m_coarse = translate(m_coarse, [0,0,height])
    # m_conf = weld(m_open, m_coarse)

    return m_open, m_top
end


"""
    meshcube(size, h) -> mesh

Mesh a conforming cube.
"""
meshcube(size, h) =  weld(mesh_box_width_lid(size,h,h)...)
