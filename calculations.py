from imports import *

def main():
    # parameter
    Aplus = 0.5
    Across = 0
    omega = 1
    phi = 0
    psi = 0
    level = 5
    N = 100
    T = 5*2*np.pi / omega

    # plotting
    graviwave = GraviWave(Aplus=Aplus, Across=Across, omega=omega, phi=phi, psi=psi)
    x = np.linspace(-5, 5, N)
    y = np.linspace(-5, 5, N)
    X, Y = np.meshgrid(x, y)
    Z = graviwave.infinitesimal_distance(t=0, x=X, y=Y, z=0)

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect="equal")
    fig.subplots_adjust(bottom=0.3)
    ax.contour(X, Y, Z, levels=[level])

    axT = plt.axes([0.25, 0.15, 0.45, 0.03])
    tSlider = Slider(axT, r"$t$", 0, T, valinit=0,valfmt="%1.1f")

    def update(val):
        t = tSlider.val
        Z = graviwave.infinitesimal_distance(t=t, x=X, y=Y, z=0)
        ax.cla()
        ax.contour(X, Y, Z, levels=[level])
        fig.canvas.draw_idle()

    tSlider.on_changed(update)

    plt.show()

if __name__ == "__main__":
    main()
