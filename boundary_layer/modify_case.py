def modify_init(u2d, v2d, k2d, om2d, vis2d):
    # set inlet field in entre domain
    u2d = np.repeat(u_bc_west[None, :], repeats=ni, axis=0)
    k2d = np.repeat(k_bc_west[None, :], repeats=ni, axis=0)
    om2d = np.repeat(om_bc_west[None, :], repeats=ni, axis=0)

    vis2d = k2d / om2d + viscos

    return u2d, v2d, k2d, om2d, vis2d


def modify_inlet():
    global y_rans, y_rans, u_rans, v_rans, k_rans, om_rans, uv_rans, k_bc_west, eps_bc_west, om_bc_west

    return u_bc_west, v_bc_west, k_bc_west, om_bc_west, u2d_face_w, convw


def modify_conv(convw, convs):
    convs[:, 0, :] = 0
    convs[:, -1, :] = 0

    return convw, convs


def modify_u(su2d, sp2d):
    global file1

    su2d[0, :] = su2d[0, :] + convw[0, :] * u_bc_west
    sp2d[0, :] = sp2d[0, :] - convw[0, :]
    vist = (
        vis2d[
            0,
            :,
        ]
        - viscos
    )
    su2d[0, :] = su2d[0, :] + vist * aw_bound * u_bc_west
    sp2d[0, :] = sp2d[0, :] - vist * aw_bound

    if iter == 0:
        print(
            "u(5)=%7.3E,u(10)=%7.3E,u(20)=%7.3E,u(30)=%7.3E,u(40)=%7.3E,u(50)=%7.3E,u(60)=%7.3E"
            % (
                u2d[ni - 5, 5],
                u2d[ni - 5, 10],
                u2d[ni - 5, 20],
                u2d[ni - 5, 30],
                u2d[ni - 5, 40],
                u2d[ni - 5, 50],
                u2d[ni - 5, 60],
            )
        )

    if iter == 0:
        print("file1 opened")
        l1 = [
            iter,
            u2d[ni - 5, 5],
            u2d[ni - 5, 10],
            u2d[ni - 5, 20],
            u2d[ni - 5, 30],
            u2d[ni - 5, 40],
            u2d[ni - 5, 50],
            u2d[ni - 5, 60],
        ]
        np.savetxt("u-time-history.dat", l1, newline=" ")
        file1 = open("u-time-history.dat", "a")  # append
    else:
        print("file1 printed")
        file1.write("\n")
        l1 = [
            iter,
            u2d[ni - 5, 5],
            u2d[ni - 5, 10],
            u2d[ni - 5, 20],
            u2d[ni - 5, 30],
            u2d[ni - 5, 40],
            u2d[ni - 5, 50],
            u2d[ni - 5, 60],
        ]
        np.savetxt(file1, l1, newline=" ")
    #         file1.write("\n")

    return su2d, sp2d


def modify_v(su2d, sp2d):
    su2d[0, :] = su2d[0, :] + convw[0, :] * v_bc_west
    sp2d[0, :] = sp2d[0, :] - convw[0, :]
    vist = vis2d[0, :] - viscos
    su2d[0, :] = su2d[0, :] + vist * aw_bound * v_bc_west
    sp2d[0, :] = sp2d[0, :] - vist * aw_bound

    return su2d, sp2d


def modify_k(su2d, sp2d):
    su2d[0, :] = su2d[0, :] + np.maximum(convw[0, :], 0) * k_bc_west
    sp2d[0, :] = sp2d[0, :] - convw[0, :]
    vist = vis2d[0, :] - viscos
    su2d[0, :] = su2d[0, :] + vist * aw_bound * k_bc_west
    sp2d[0, :] = sp2d[0, :] - vist * aw_bound

    return su2d, sp2d


def modify_om(su2d, sp2d):
    su2d[0, :] = su2d[0, :] + np.maximum(convw[0, :], 0) * om_bc_west
    sp2d[0, :] = sp2d[0, :] - convw[0, :]
    vist = vis2d[0, :] - viscos
    su2d[0, :] = su2d[0, :] + vist * aw_bound * om_bc_west
    sp2d[0, :] = sp2d[0, :] - vist * aw_bound

    return su2d, sp2d


def modify_outlet(convw):
    # inlet
    flow_in = np.sum(convw[0, :])
    flow_out = np.sum(convw[-1, :])
    #  flow_out=np.sum(convw[-2,:])
    area_out = np.sum(areaw[-1, :])

    uinc = (flow_in - flow_out) / area_out
    ares = areaw[-1, :]
    convw[-1, :] = convw[-1, :] + uinc * ares
    #  convw[-1,:]=convw[-2,:]+uinc*ares

    print("area_out", area_out)

    flow_out_new = np.sum(convw[-1, :])

    print(
        f"{'flow_in: '} {flow_in:.3e},{'  flow_out: '} {flow_out:.3e},{'  flow_out_new: '} {flow_out_new:.3e},{'  uinc: '} {uinc:.3e}"
    )

    return convw


def fix_omega():
    aw2d[:, 0] = 0
    ae2d[:, 0] = 0
    as2d[:, 0] = 0
    an2d[:, 0] = 0
    ap2d[:, 0] = 1
    su2d[:, 0] = om_bc_south

    return aw2d, ae2d, as2d, an2d, ap2d, su2d, sp2d


def modify_vis(vis2d):
    return vis2d


def fix_k():
    return aw2d, ae2d, as2d, an2d, ap2d, su2d, sp2d
