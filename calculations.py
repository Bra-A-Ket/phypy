from imports import *

def main():
    graviwave = GraviWave(Aplus=0.5, Across=0)
    x = np.linspace(-5, 5, 100)
    y = np.linspace(-5, 5, 100)
    X, Y = np.meshgrid(x, y)
    Z = graviwave.infinitesimal_distance(t=0, x=X, y=Y, z=0)


    fig = plt.figure()
    ax = fig.add_subplot(111, aspect="equal")
    fig.subplots_adjust(bottom=0.3)
    ax.contour(X, Y, Z, levels=[5])

    axT = plt.axes([0.25, 0.15, 0.45, 0.03])
    tSlider = Slider(axT, r"$t$", 0, 15, valinit=0,valfmt="%1.1f")

    def update(val):
        t = tSlider.val
        Z = graviwave.infinitesimal_distance(t=t, x=X, y=Y, z=0)
        ax.cla()
        ax.contour(X, Y, Z, levels=[5])
        fig.canvas.draw_idle()

    tSlider.on_changed(update)

    plt.show()

if __name__ == "__main__":
    main()
