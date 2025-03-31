from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
from OCC.Core.BRepTools import breptools
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TopoDS import topods
from OCC.Core.gp import gp_Pnt
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_VERTEX, TopAbs_REVERSED
from OCC.Core.BRep import BRep_Tool
from OCC.Extend.TopologyUtils import TopologyExplorer


def create_box_from_corners(xmin, ymin, zmin, xmax, ymax, zmax):
    # 计算长宽高
    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin

    # 创建原点（最小角点）
    origin = gp_Pnt(xmin, ymin, zmin)

    # 生成长方体
    box = BRepPrimAPI_MakeBox(origin, dx, dy, dz).Shape()

    return box

# 示例用法
if __name__ == "__main__":
    # 输入对角点坐标（支持任意顺序，函数内部自动处理）
    box_shape = create_box_from_corners(500, 151, -2187.5, 6900, 2165, 37.5)
    topo_exp = TopologyExplorer(box_shape)
    off = 1900759

    # Lists to store vertices, normals, and faces for each group
    all_vertices = []
    all_normals = []
    face_meshes = []
    vertex_offset = 0
    normal_offset = 0

    for i, face in enumerate(topo_exp.faces()):

        breptools.Clean(face)
        # Discretize face to a mesh
        BRepMesh_IncrementalMesh(face, 0.2, True, 0.36).Perform()  # Ensure meshing for each face
        rev_flag = False
        if face.Orientation() == TopAbs_REVERSED:
            rev_flag = True

        # Get the triangulation for this face
        loc = TopLoc_Location()  # Initialize a TopLoc_Location object
        face_triangulation = BRep_Tool.Triangulation(topods.Face(face), loc)

        if face_triangulation is None:
            continue  # Skip faces without triangulation

        # Extract vertices and faces from the triangulation
        vertices = []
        faces = []

        for j in range(1, face_triangulation.NbNodes() + 1):
            vertex = face_triangulation.Node(j)
            vertices.append([vertex.X(), vertex.Y(), vertex.Z()])

        for j in range(1, face_triangulation.NbTriangles() + 1):
            triangle = face_triangulation.Triangle(j)
            # Get the vertices for the triangle
            if rev_flag:
                v1_idx = triangle.Value(3) - 1
                v2_idx = triangle.Value(2) - 1
                v3_idx = triangle.Value(1) - 1
            else:
                v1_idx = triangle.Value(1) - 1
                v2_idx = triangle.Value(2) - 1
                v3_idx = triangle.Value(3) - 1
            v1 = vertices[v1_idx]
            v2 = vertices[v2_idx]
            v3 = vertices[v3_idx]

            # Compute the normal for this triangle



            # Store face with vertex indices and corresponding normal
            face_vertex_indices = [v1_idx + vertex_offset, v2_idx + vertex_offset, v3_idx + vertex_offset]
            faces.append(face_vertex_indices)

        # Store the face mesh, vertices, and normals
        face_meshes.append((f"face{i}", faces))
        all_vertices.extend(vertices)
        vertex_offset += len(vertices)

    # Write to OBJ file with vertices and normals first, then face groups
    with open("outbox.obj", 'w') as obj:
        # Write all vertices at the beginning
        for v in all_vertices:
            obj.write(f"v {v[0]} {v[1]} {v[2] }\n")

        obj.write("\n")


        obj.write("\n")

        # Write each face group
        for group_name, faces in face_meshes:
            obj.write(f"g box{group_name}\n")

            for f_verts in faces:
                obj.write(
                    f"f {f_verts[0] + off}// {f_verts[1] + off}// {f_verts[2] + off}//\n")
            obj.write("\n")

    print(f"Conversion complete: outbox.ob")