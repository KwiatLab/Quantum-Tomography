
def displayOutput(p,inten,fval,errs,mean,mean_fid):
    for x in range(0, len(err_functions)):
        print(err_functions[x], "=", mean[x], "+/-", errs[x])
    print('Mean fidelity =', mean_fid)
    print("Intensity =", inten)
    print("fval =", fval)
    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    xpos = [.5, .5, .5, .5,
            1.5, 1.5, 1.5, 1.5,
            2.5, 2.5, 2.5, 2.5,
            3.5, 3.5, 3.5, 3.5]
    ypos = [.5, 1.5, 2.5, 3.5,
            .5, 1.5, 2.5, 3.5,
            .5, 1.5, 2.5, 3.5,
            .5, 1.5, 2.5, 3.5, ]
    zpos = np.zeros(16)
    dx = [.9, .9, .9, .9
        , .9, .9, .9, .9
        , .9, .9, .9, .9
        , .9, .9, .9, .9]
    dy = [.9, .9, .9, .9
        , .9, .9, .9, .9
        , .9, .9, .9, .9
        , .9, .9, .9, .9]

    dz = [.5, 1.5, 2.5, 3.5,
          .5, 1.5, 2.5, 3.5,
          .5, 1.5, 2.5, 3.5,
          .5, 1.5, 2.5, 3.5, ]
    dz = p.flatten().astype(float)


    colors = ['r','g','b','purple']

    #ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, color='#00ceaa')
    for i in range(4):
        ax1.bar3d(xpos[i*4:i*4+4], ypos[i*4:i*4+4], zpos[i*4:i*4+4], dx[i*4:i*4+4], dy[i*4:i*4+4], dz[i*4:i*4+4],edgecolor="black",alpha=.65, color=colors[i])

    ax1.axes.set_xticklabels(["|HH>", "|HV>", "|VH>", "|VV>"])
    ax1.axes.set_yticklabels(["<HH|", "<HV|", "<VH|", "<VV|"])
    ax1.axes.set_xticks([1, 2, 3, 4])
    ax1.axes.set_yticks([1, 2, 3, 4])
    ax1.axes.set_zticks([-1.0, -.8, -.6, -.4, -.2, 0, .2, .4, .6, .8, 1.0])
    ax1.axes.set_zlim3d(-1, 1)
    plt.title("Rho")



    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    xpos = [.5, .5, .5, .5,
            1.5, 1.5, 1.5, 1.5,
            2.5, 2.5, 2.5, 2.5,
            3.5, 3.5, 3.5, 3.5]
    ypos = [.5, 1.5, 2.5, 3.5,
            .5, 1.5, 2.5, 3.5,
            .5, 1.5, 2.5, 3.5,
            .5, 1.5, 2.5, 3.5, ]
    zpos = np.zeros(16)
    dx = [.9, .9, .9, .9
        , .9, .9, .9, .9
        , .9, .9, .9, .9
        , .9, .9, .9, .9]
    dy = [.9, .9, .9, .9
        , .9, .9, .9, .9
        , .9, .9, .9, .9
        , .9, .9, .9, .9]

    dz = [.5, 1.5, 2.5, 3.5,
          .5, 1.5, 2.5, 3.5,
          .5, 1.5, 2.5, 3.5,
          .5, 1.5, 2.5, 3.5, ]
    dz = p.flatten().imag.astype(float)
    colors = ['r', 'g', 'b', 'purple']

    # ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, color='#00ceaa')
    for i in range(4):
        ax1.bar3d(xpos[i * 4:i * 4 + 4], ypos[i * 4:i * 4 + 4], zpos[i * 4:i * 4 + 4], dx[i * 4:i * 4 + 4],
                  dy[i * 4:i * 4 + 4], dz[i * 4:i * 4 + 4], edgecolor="black", alpha=.65, color=colors[i])


    ax1.axes.set_xticklabels(["|HH>", "|HV>", "|VH>", "|VV>"])
    ax1.axes.set_yticklabels(["<HH|", "<HV|", "<VH|", "<VV|"])
    ax1.axes.set_xticks([1, 2, 3, 4])
    ax1.axes.set_yticks([1, 2, 3, 4])
    ax1.axes.set_zticks([-1.0, -.8, -.6, -.4, -.2, 0, .2, .4, .6, .8, 1.0])
    ax1.axes.set_zlim3d(-1, 1)
    plt.title("Rho Imaginary")
    plt.show()
