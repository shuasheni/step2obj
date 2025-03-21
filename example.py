from step2mesh.box_processer import BounderBoxProcessor


if __name__ == "__main__":
    # 4083028面 => 2200364面

    processor = BounderBoxProcessor(step_file="501000000797UTF.STEP", show_detail=True, scale=0.9, p=0.1)
    assert processor.travel(), "not assembly"
    processor.cal_cutbox()
    print(f"cut box: ({processor.cut_box[0]}, {processor.cut_box[1]}, {processor.cut_box[2]}), "
          f"({processor.cut_box[3]}, {processor.cut_box[4]}, {processor.cut_box[5]})")
    all_vertices = [
        [processor.cut_box[0], processor.cut_box[1], processor.cut_box[2]],
        [processor.cut_box[0], processor.cut_box[1], processor.cut_box[5]],
        [processor.cut_box[0], processor.cut_box[4], processor.cut_box[2]],
        [processor.cut_box[0], processor.cut_box[4], processor.cut_box[5]],
        [processor.cut_box[3], processor.cut_box[1], processor.cut_box[2]],
        [processor.cut_box[3], processor.cut_box[1], processor.cut_box[5]],
        [processor.cut_box[3], processor.cut_box[4], processor.cut_box[2]],
        [processor.cut_box[3], processor.cut_box[4], processor.cut_box[5]],
    ]
    all_face_meshes = [
        # 底面
        [[0, 1, 3], [0, 3, 2]],
        # 顶面
        [[4, 5, 7], [4, 7, 6]],
        # 前面
        [[0, 1, 5], [0, 5, 4]],
        # 后面
        [[2, 3, 7], [2, 7, 6]],
        # 左面
        [[0, 2, 6], [0, 6, 4]],
        # 右面
        [[1, 3, 7], [1, 7, 5]]
    ]
    box_info = ('box', None, all_vertices, all_face_meshes, processor.cut_box)

    processor.final_shapes.append(box_info)

    processor.compute_convex_hull_all()
    fnum = processor.output(obj_filename="test123432123.obj")






