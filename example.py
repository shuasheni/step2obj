from step2mesh.box_processer import BounderBoxProcessor


if __name__ == "__main__":
    # 4083028面 => 2200364面

    processor = BounderBoxProcessor(step_file="501000000797UTF.STEP", show_detail=True, scale=0.8)
    assert processor.travel(), "not assembly"

    print(f"bounding box: ({processor.bounder_box[0]}, {processor.bounder_box[1]}, {processor.bounder_box[2]}), "
          f"({processor.bounder_box[3]}, {processor.bounder_box[4]}, {processor.bounder_box[5]})")
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
        [[1, 2, 4], [1, 4, 3]],
        # 顶面
        [[5, 6, 8], [5, 8, 7]],
        # 前面
        [[1, 2, 6], [1, 6, 5]],
        # 后面
        [[3, 4, 8], [3, 8, 7]],
        # 左面
        [[1, 3, 7], [1, 7, 5]],
        # 右面
        [[2, 4, 8], [2, 8, 6]]
    ]
    box_info = ('box', None, all_vertices, all_face_meshes, processor.cut_box)

    processor.final_shapes.append(box_info)

    processor.compute_convex_hull_all()
    fnum = processor.output(obj_filename="test12343212.obj")






