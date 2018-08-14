#include <QMainWindow>
#include <QImage>
#include <QCloseEvent>
#include <QResizeEvent>
#include <QMouseEvent>
#include <QKeyEvent>
#include <ctime>
#include <vector>
#include "ui_mainwindow.h"
#include "testwindow.h"
#include "helperfunctions.h"

namespace Ui {
class MainWindow;
}

namespace smallpaint_fixed {
struct Vec;
}

namespace smallpaint_painterly {
struct Vec;
}

namespace smallpaint_bvh {
struct Vec;
}

namespace smallpaint_pssmlt {
struct Vec;
}

namespace smallpaint_smallmedia {
struct Vec;
}

namespace smallpaint_ppm {
struct Vec;
}

namespace smallpaint_vrl {
struct Vec;
}

class MainWindow : public QMainWindow {
	Q_OBJECT

		friend class HelperFunctions;

public:
	int width, height;
	explicit MainWindow(QWidget *parent = 0);
	~MainWindow();
	void drawImage(QImage img, std::string name, int currentSpp, int goalSpp);

public slots:
	void on_renderButton_clicked();
    void saveImage();

protected:
	void closeEvent(QCloseEvent *event);
	void resizeEvent(QResizeEvent *event);
	void keyPressEvent(QKeyEvent *event);

signals:
	void smallpaint_fixedSignal(smallpaint_fixed::Vec **pix, int spp);
	void smallpaint_painterlySignal(smallpaint_painterly::Vec **pix, int spp);
	void smallpaint_bvhSignal(smallpaint_bvh::Vec **pix, int spp);
	void smallpaint_pssmltSignal(smallpaint_pssmlt::Vec **pix, int spp);
	void smallpaint_smallmediaSignal(smallpaint_smallmedia::Vec **pix, int spp);
    void smallpaint_ppmSignal(smallpaint_ppm::Vec **pix, int spp);
    void smallpaint_vrlSignal(smallpaint_vrl::Vec **pix, int spp);

	private slots:
	void smallpaint_fixedSlot(smallpaint_fixed::Vec **pix, int spp);
	void smallpaint_painterlySlot(smallpaint_painterly::Vec **pix, int spp);
	void smallpaint_bvhSlot(smallpaint_bvh::Vec **pix, int spp);
	void smallpaint_pssmltSlot(smallpaint_pssmlt::Vec **pix, int spp);
	void smallpaint_smallmediaSlot(smallpaint_smallmedia::Vec **pix, int spp);
    void smallpaint_ppmSlot(smallpaint_ppm::Vec **pix, int spp);
    void smallpaint_vrlSlot(smallpaint_vrl::Vec **pix, int spp);
    void on_renderModesComboBox_currentTextChanged(const QString &arg1);

    void on_vrlScene_currentIndexChanged(int index);

private:
	Ui::MainWindow *ui;
	HelperFunctions helper;

    QImage image = QImage("../test_images/logo.png");
	std::string lastImage = "";
	const int inf = 1e9;
	clock_t start;

	TestWindow *testWindow;
	std::vector<int> tests;
	bool testing = false;

	void resizeImage();
	void testRender(bool renderButtonActivated);
};

void setInstanceForRenderDistribution(MainWindow* w);
void setInstanceForHelperFunctions(MainWindow* w);
