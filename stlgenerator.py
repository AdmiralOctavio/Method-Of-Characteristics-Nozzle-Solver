import numpy as np
from stl import mesh

def create_stl(x_contour, y_contour, Mach):
    filename = f"Nozzle_Contour_M={Mach:.2f}.stl"
    num_points = len(x_contour)
    num_facets = 60
    theta = np.linspace(0, 2 * np.pi, num_facets, endpoint=False) 

    num_triangles = 2 * (num_points - 1) * num_facets
    
    # Create the mesh array to hold vertices (3 points per triangle, 3 coordinates per point)
    nozzle_mesh = np.zeros(num_triangles, dtype=mesh.Mesh.dtype)

    # 2. GENERATE VERTICES (Revolution)
    # The vertices are arranged in rings along the contour
    vertices = np.zeros((num_points, num_facets, 3))
    
    for i in range(num_points):
        # x is constant for each ring (axial position)
        x_val = x_contour[i]
        # y is cosine of theta * radius
        # z is sine of theta * radius
        r_val = y_contour[i]
        
        # Calculate vertices for the current ring (x, y, z)
        vertices[i, :, 0] = x_val
        vertices[i, :, 1] = r_val * np.cos(theta)
        vertices[i, :, 2] = r_val * np.sin(theta)

    # 3. GENERATE TRIANGLES (Faces)
    triangle_index = 0
    
    # Iterate through the segments along the nozzle length (i)
    for i in range(num_points - 1):
        # Iterate through the facets around the circumference (j)
        for j in range(num_facets):
            # Define vertices for the current quad (two triangles)
            # Current ring (i)
            v1 = vertices[i, j]
            v2 = vertices[i, (j + 1) % num_facets] # Wrap around to the first facet
            
            # Next ring (i+1)
            v3 = vertices[i + 1, (j + 1) % num_facets]
            v4 = vertices[i + 1, j]

            # Triangle 1: v1, v2, v3 (Front/leading triangle)
            nozzle_mesh['vectors'][triangle_index] = [v1, v2, v3]
            triangle_index += 1
            
            # Triangle 2: v1, v3, v4 (Rear/lagging triangle)
            nozzle_mesh['vectors'][triangle_index] = [v1, v3, v4]
            triangle_index += 1

    # 4. SAVE MESH
    # Create the mesh object and save it to the file
    mesh_obj = mesh.Mesh(nozzle_mesh)
    
    # You might want to flip the normals if the surface appears inverted in CAD:
    # mesh_obj.flip_normals() 
    
    mesh_obj.save(filename)

