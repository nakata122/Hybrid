/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.7.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QFormLayout>
#include <QtWidgets/QFrame>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QScrollArea>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QStackedWidget>
#include <QtWidgets/QTextBrowser>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include "image.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralWidget;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout;
    QFrame *frame;
    QVBoxLayout *modes;
    QComboBox *renderModesComboBox;
    QStackedWidget *renderModesSettings;
    QWidget *fixed;
    QVBoxLayout *verticalLayout_3;
    QHBoxLayout *horizontalLayout_2;
    QLabel *fixedSizeLabel;
    QLineEdit *fixedSize;
    QLabel *label;
    QFrame *line;
    QHBoxLayout *horizontalLayout_3;
    QLabel *fixedSamplesLabel;
    QLineEdit *fixedSamples;
    QLabel *label_2;
    QFrame *line_2;
    QHBoxLayout *horizontalLayout_4;
    QLabel *fixedRefractionLabel;
    QLineEdit *fixedRefraction;
    QLabel *label_3;
    QSpacerItem *verticalSpacer;
    QWidget *painterly;
    QVBoxLayout *verticalLayout_4;
    QHBoxLayout *horizontalLayout_7;
    QLabel *painterlySizeLabel;
    QLineEdit *painterlySize;
    QLabel *label_4;
    QFrame *line_3;
    QHBoxLayout *horizontalLayout_5;
    QLabel *painterlySamplesLabel;
    QLineEdit *painterlySamples;
    QLabel *label_5;
    QFrame *line_4;
    QHBoxLayout *horizontalLayout_6;
    QLabel *painterlyRefractionLabel;
    QLineEdit *painterlyRefraction;
    QLabel *label_6;
    QSpacerItem *verticalSpacer_2;
    QWidget *bvh;
    QVBoxLayout *verticalLayout_5;
    QHBoxLayout *horizontalLayout_10;
    QLabel *bvhSizeLabel;
    QLineEdit *bvhSize;
    QLabel *label_7;
    QFrame *line_5;
    QHBoxLayout *horizontalLayout_9;
    QLabel *bvhSamplesLabel;
    QLineEdit *bvhSamples;
    QLabel *label_8;
    QFrame *line_6;
    QHBoxLayout *horizontalLayout_8;
    QLabel *bvhRefractionLabel;
    QLineEdit *bvhRefraction;
    QLabel *label_9;
    QSpacerItem *verticalSpacer_3;
    QWidget *pssmlt;
    QVBoxLayout *verticalLayout_6;
    QHBoxLayout *horizontalLayout_11;
    QLabel *pssmltSizeLabel;
    QLineEdit *pssmltSize;
    QLabel *label_10;
    QFrame *line_7;
    QHBoxLayout *horizontalLayout_12;
    QLabel *pssmltSamplesLabel;
    QLineEdit *pssmltSamples;
    QLabel *label_11;
    QFrame *line_8;
    QHBoxLayout *horizontalLayout_13;
    QLabel *pssmltRefractionLabel;
    QLineEdit *pssmltRefraction;
    QLabel *label_12;
    QSpacerItem *verticalSpacer_4;
    QWidget *smallmedia;
    QVBoxLayout *verticalLayout_7;
    QHBoxLayout *horizontalLayout_14;
    QLabel *smallmediaSizeLabel;
    QLineEdit *smallmediaSize;
    QLabel *label_13;
    QFrame *line_9;
    QHBoxLayout *horizontalLayout_15;
    QLabel *smallmediaSamplesLabel;
    QLineEdit *smallmediaSamples;
    QLabel *label_14;
    QFrame *line_10;
    QHBoxLayout *horizontalLayout_16;
    QLabel *smallmediaRefractionLabel;
    QLineEdit *smallmediaRefraction;
    QLabel *label_15;
    QSpacerItem *verticalSpacer_5;
    QWidget *ppm;
    QVBoxLayout *verticalLayout_8;
    QHBoxLayout *horizontalLayout_17;
    QLabel *ppmSizeLabel;
    QLineEdit *ppmSize;
    QLabel *label_16;
    QFrame *line_11;
    QHBoxLayout *horizontalLayout_18;
    QLabel *ppmSamplesLabel;
    QLineEdit *ppmSamples;
    QLabel *label_17;
    QFrame *line_12;
    QHBoxLayout *horizontalLayout_19;
    QLabel *ppmRefractionLabel;
    QLineEdit *ppmRefraction;
    QLabel *label_18;
    QSpacerItem *verticalSpacer_6;
    QWidget *vrl;
    QVBoxLayout *verticalLayout_2;
    QScrollArea *scrollArea;
    QWidget *scrollAreaWidgetContents;
    QVBoxLayout *verticalLayout_9;
    QHBoxLayout *horizontalLayout_22;
    QLabel *vrlSizeLabel;
    QLineEdit *vrlSize;
    QLabel *label_19;
    QFrame *line_13;
    QHBoxLayout *horizontalLayout_20;
    QLabel *vrlSamplesLabel;
    QLineEdit *vrlSamples;
    QLabel *label_20;
    QFrame *line_14;
    QHBoxLayout *horizontalLayout_23;
    QLabel *vrlVRLPSLabel;
    QLineEdit *vrlVRLPS;
    QLabel *label_23;
    QFrame *line_16;
    QHBoxLayout *horizontalLayout_26;
    QLabel *vrlBouncesLabel;
    QLineEdit *vrlBounces;
    QLabel *label_27;
    QFrame *line_17;
    QHBoxLayout *horizontalLayout_27;
    QLabel *vrlRefractionLabel;
    QLineEdit *vrlRefraction;
    QLabel *label_29;
    QFrame *line_19;
    QHBoxLayout *horizontalLayout_25;
    QLabel *vrlMediumLabel;
    QLineEdit *vrlSigmaA;
    QLineEdit *vrlSigmaS;
    QLineEdit *vrlG;
    QLabel *label_25;
    QFrame *line_21;
    QHBoxLayout *horizontalLayout_28;
    QLabel *vrlIntensityLabel;
    QLineEdit *vrlIntensity;
    QLabel *label_26;
    QSpacerItem *verticalSpacer_7;
    QFrame *line_20;
    QCheckBox *vrlMediumRadiance;
    QCheckBox *vrlSurfaceRadiance;
    QFrame *line_15;
    QHBoxLayout *horizontalLayout_21;
    QLabel *label_21;
    QComboBox *vrlSampling;
    QLabel *label_22;
    QFrame *line_18;
    QHBoxLayout *horizontalLayout_24;
    QLabel *label_24;
    QComboBox *vrlScene;
    QWidget *tests;
    QFormLayout *formLayout_5;
    QRadioButton *fixedTest;
    QRadioButton *painterlyTest;
    QRadioButton *bvhTest;
    QRadioButton *pssmltTest;
    QRadioButton *smallmediaTest;
    QRadioButton *ppmTest;
    QRadioButton *vrlTest;
    QTextBrowser *timeInfo;
    QProgressBar *progressBar;
    Image *image;
    QHBoxLayout *horizontalLayoutBottom;
    QTextBrowser *renderInfo;
    QPushButton *renderButton;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QStringLiteral("MainWindow"));
        MainWindow->resize(800, 579);
        QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(MainWindow->sizePolicy().hasHeightForWidth());
        MainWindow->setSizePolicy(sizePolicy);
        MainWindow->setMinimumSize(QSize(800, 579));
        MainWindow->setMaximumSize(QSize(16777215, 16777215));
        MainWindow->setStyleSheet(QStringLiteral(""));
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        centralWidget->setMinimumSize(QSize(800, 579));
        verticalLayout = new QVBoxLayout(centralWidget);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setSpacing(6);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        horizontalLayout->setSizeConstraint(QLayout::SetMaximumSize);
        frame = new QFrame(centralWidget);
        frame->setObjectName(QStringLiteral("frame"));
        frame->setMaximumSize(QSize(300, 16777215));
        frame->setSizeIncrement(QSize(0, 0));
        frame->setAutoFillBackground(true);
        frame->setStyleSheet(QStringLiteral(""));
        frame->setFrameShape(QFrame::Box);
        frame->setFrameShadow(QFrame::Raised);
        modes = new QVBoxLayout(frame);
        modes->setSpacing(3);
        modes->setContentsMargins(11, 11, 11, 11);
        modes->setObjectName(QStringLiteral("modes"));
        modes->setSizeConstraint(QLayout::SetDefaultConstraint);
        modes->setContentsMargins(3, 3, 3, 3);
        renderModesComboBox = new QComboBox(frame);
        renderModesComboBox->setObjectName(QStringLiteral("renderModesComboBox"));
        renderModesComboBox->setMaximumSize(QSize(300, 16777215));
        QFont font;
        font.setPointSize(10);
        font.setBold(false);
        font.setWeight(50);
        font.setKerning(true);
        renderModesComboBox->setFont(font);
        renderModesComboBox->setMaxVisibleItems(10);

        modes->addWidget(renderModesComboBox);

        renderModesSettings = new QStackedWidget(frame);
        renderModesSettings->setObjectName(QStringLiteral("renderModesSettings"));
        QSizePolicy sizePolicy1(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(renderModesSettings->sizePolicy().hasHeightForWidth());
        renderModesSettings->setSizePolicy(sizePolicy1);
        renderModesSettings->setMaximumSize(QSize(300, 16777215));
        QFont font1;
        font1.setKerning(true);
        renderModesSettings->setFont(font1);
        renderModesSettings->setAutoFillBackground(true);
        renderModesSettings->setFrameShape(QFrame::NoFrame);
        renderModesSettings->setFrameShadow(QFrame::Plain);
        renderModesSettings->setLineWidth(0);
        fixed = new QWidget();
        fixed->setObjectName(QStringLiteral("fixed"));
        QSizePolicy sizePolicy2(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(fixed->sizePolicy().hasHeightForWidth());
        fixed->setSizePolicy(sizePolicy2);
        fixed->setMaximumSize(QSize(16777215, 16777215));
        verticalLayout_3 = new QVBoxLayout(fixed);
        verticalLayout_3->setSpacing(6);
        verticalLayout_3->setContentsMargins(11, 11, 11, 11);
        verticalLayout_3->setObjectName(QStringLiteral("verticalLayout_3"));
        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        fixedSizeLabel = new QLabel(fixed);
        fixedSizeLabel->setObjectName(QStringLiteral("fixedSizeLabel"));
        QFont font2;
        font2.setBold(true);
        font2.setWeight(75);
        fixedSizeLabel->setFont(font2);
        fixedSizeLabel->setToolTipDuration(-1);

        horizontalLayout_2->addWidget(fixedSizeLabel);

        fixedSize = new QLineEdit(fixed);
        fixedSize->setObjectName(QStringLiteral("fixedSize"));
        fixedSize->setMaximumSize(QSize(100, 20));

        horizontalLayout_2->addWidget(fixedSize);


        verticalLayout_3->addLayout(horizontalLayout_2);

        label = new QLabel(fixed);
        label->setObjectName(QStringLiteral("label"));
        label->setMaximumSize(QSize(16777215, 16777215));
        label->setWordWrap(true);

        verticalLayout_3->addWidget(label);

        line = new QFrame(fixed);
        line->setObjectName(QStringLiteral("line"));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);

        verticalLayout_3->addWidget(line);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        fixedSamplesLabel = new QLabel(fixed);
        fixedSamplesLabel->setObjectName(QStringLiteral("fixedSamplesLabel"));
        fixedSamplesLabel->setFont(font2);

        horizontalLayout_3->addWidget(fixedSamplesLabel);

        fixedSamples = new QLineEdit(fixed);
        fixedSamples->setObjectName(QStringLiteral("fixedSamples"));
        fixedSamples->setMaximumSize(QSize(100, 20));

        horizontalLayout_3->addWidget(fixedSamples);


        verticalLayout_3->addLayout(horizontalLayout_3);

        label_2 = new QLabel(fixed);
        label_2->setObjectName(QStringLiteral("label_2"));
        label_2->setWordWrap(true);

        verticalLayout_3->addWidget(label_2);

        line_2 = new QFrame(fixed);
        line_2->setObjectName(QStringLiteral("line_2"));
        line_2->setFrameShape(QFrame::HLine);
        line_2->setFrameShadow(QFrame::Sunken);

        verticalLayout_3->addWidget(line_2);

        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setSpacing(6);
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        fixedRefractionLabel = new QLabel(fixed);
        fixedRefractionLabel->setObjectName(QStringLiteral("fixedRefractionLabel"));
        fixedRefractionLabel->setFont(font2);

        horizontalLayout_4->addWidget(fixedRefractionLabel);

        fixedRefraction = new QLineEdit(fixed);
        fixedRefraction->setObjectName(QStringLiteral("fixedRefraction"));
        fixedRefraction->setMaximumSize(QSize(100, 20));

        horizontalLayout_4->addWidget(fixedRefraction);


        verticalLayout_3->addLayout(horizontalLayout_4);

        label_3 = new QLabel(fixed);
        label_3->setObjectName(QStringLiteral("label_3"));
        label_3->setWordWrap(true);

        verticalLayout_3->addWidget(label_3);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_3->addItem(verticalSpacer);

        renderModesSettings->addWidget(fixed);
        painterly = new QWidget();
        painterly->setObjectName(QStringLiteral("painterly"));
        verticalLayout_4 = new QVBoxLayout(painterly);
        verticalLayout_4->setSpacing(6);
        verticalLayout_4->setContentsMargins(11, 11, 11, 11);
        verticalLayout_4->setObjectName(QStringLiteral("verticalLayout_4"));
        horizontalLayout_7 = new QHBoxLayout();
        horizontalLayout_7->setSpacing(6);
        horizontalLayout_7->setObjectName(QStringLiteral("horizontalLayout_7"));
        painterlySizeLabel = new QLabel(painterly);
        painterlySizeLabel->setObjectName(QStringLiteral("painterlySizeLabel"));
        painterlySizeLabel->setFont(font2);

        horizontalLayout_7->addWidget(painterlySizeLabel);

        painterlySize = new QLineEdit(painterly);
        painterlySize->setObjectName(QStringLiteral("painterlySize"));
        painterlySize->setMaximumSize(QSize(100, 20));

        horizontalLayout_7->addWidget(painterlySize);


        verticalLayout_4->addLayout(horizontalLayout_7);

        label_4 = new QLabel(painterly);
        label_4->setObjectName(QStringLiteral("label_4"));
        label_4->setWordWrap(true);

        verticalLayout_4->addWidget(label_4);

        line_3 = new QFrame(painterly);
        line_3->setObjectName(QStringLiteral("line_3"));
        line_3->setFrameShape(QFrame::HLine);
        line_3->setFrameShadow(QFrame::Sunken);

        verticalLayout_4->addWidget(line_3);

        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setSpacing(6);
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        painterlySamplesLabel = new QLabel(painterly);
        painterlySamplesLabel->setObjectName(QStringLiteral("painterlySamplesLabel"));
        painterlySamplesLabel->setFont(font2);

        horizontalLayout_5->addWidget(painterlySamplesLabel);

        painterlySamples = new QLineEdit(painterly);
        painterlySamples->setObjectName(QStringLiteral("painterlySamples"));
        painterlySamples->setMaximumSize(QSize(100, 20));

        horizontalLayout_5->addWidget(painterlySamples);


        verticalLayout_4->addLayout(horizontalLayout_5);

        label_5 = new QLabel(painterly);
        label_5->setObjectName(QStringLiteral("label_5"));
        label_5->setWordWrap(true);

        verticalLayout_4->addWidget(label_5);

        line_4 = new QFrame(painterly);
        line_4->setObjectName(QStringLiteral("line_4"));
        line_4->setFrameShape(QFrame::HLine);
        line_4->setFrameShadow(QFrame::Sunken);

        verticalLayout_4->addWidget(line_4);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setSpacing(6);
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        painterlyRefractionLabel = new QLabel(painterly);
        painterlyRefractionLabel->setObjectName(QStringLiteral("painterlyRefractionLabel"));
        painterlyRefractionLabel->setFont(font2);

        horizontalLayout_6->addWidget(painterlyRefractionLabel);

        painterlyRefraction = new QLineEdit(painterly);
        painterlyRefraction->setObjectName(QStringLiteral("painterlyRefraction"));
        painterlyRefraction->setMaximumSize(QSize(100, 20));

        horizontalLayout_6->addWidget(painterlyRefraction);


        verticalLayout_4->addLayout(horizontalLayout_6);

        label_6 = new QLabel(painterly);
        label_6->setObjectName(QStringLiteral("label_6"));
        label_6->setWordWrap(true);

        verticalLayout_4->addWidget(label_6);

        verticalSpacer_2 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_4->addItem(verticalSpacer_2);

        renderModesSettings->addWidget(painterly);
        bvh = new QWidget();
        bvh->setObjectName(QStringLiteral("bvh"));
        verticalLayout_5 = new QVBoxLayout(bvh);
        verticalLayout_5->setSpacing(6);
        verticalLayout_5->setContentsMargins(11, 11, 11, 11);
        verticalLayout_5->setObjectName(QStringLiteral("verticalLayout_5"));
        horizontalLayout_10 = new QHBoxLayout();
        horizontalLayout_10->setSpacing(6);
        horizontalLayout_10->setObjectName(QStringLiteral("horizontalLayout_10"));
        bvhSizeLabel = new QLabel(bvh);
        bvhSizeLabel->setObjectName(QStringLiteral("bvhSizeLabel"));
        bvhSizeLabel->setFont(font2);

        horizontalLayout_10->addWidget(bvhSizeLabel);

        bvhSize = new QLineEdit(bvh);
        bvhSize->setObjectName(QStringLiteral("bvhSize"));
        bvhSize->setMaximumSize(QSize(100, 20));

        horizontalLayout_10->addWidget(bvhSize);


        verticalLayout_5->addLayout(horizontalLayout_10);

        label_7 = new QLabel(bvh);
        label_7->setObjectName(QStringLiteral("label_7"));
        QFont font3;
        font3.setItalic(false);
        label_7->setFont(font3);
        label_7->setWordWrap(true);

        verticalLayout_5->addWidget(label_7);

        line_5 = new QFrame(bvh);
        line_5->setObjectName(QStringLiteral("line_5"));
        line_5->setFrameShape(QFrame::HLine);
        line_5->setFrameShadow(QFrame::Sunken);

        verticalLayout_5->addWidget(line_5);

        horizontalLayout_9 = new QHBoxLayout();
        horizontalLayout_9->setSpacing(6);
        horizontalLayout_9->setObjectName(QStringLiteral("horizontalLayout_9"));
        bvhSamplesLabel = new QLabel(bvh);
        bvhSamplesLabel->setObjectName(QStringLiteral("bvhSamplesLabel"));
        bvhSamplesLabel->setFont(font2);

        horizontalLayout_9->addWidget(bvhSamplesLabel);

        bvhSamples = new QLineEdit(bvh);
        bvhSamples->setObjectName(QStringLiteral("bvhSamples"));
        bvhSamples->setMaximumSize(QSize(100, 20));

        horizontalLayout_9->addWidget(bvhSamples);


        verticalLayout_5->addLayout(horizontalLayout_9);

        label_8 = new QLabel(bvh);
        label_8->setObjectName(QStringLiteral("label_8"));
        label_8->setFont(font3);
        label_8->setWordWrap(true);

        verticalLayout_5->addWidget(label_8);

        line_6 = new QFrame(bvh);
        line_6->setObjectName(QStringLiteral("line_6"));
        line_6->setFrameShape(QFrame::HLine);
        line_6->setFrameShadow(QFrame::Sunken);

        verticalLayout_5->addWidget(line_6);

        horizontalLayout_8 = new QHBoxLayout();
        horizontalLayout_8->setSpacing(6);
        horizontalLayout_8->setObjectName(QStringLiteral("horizontalLayout_8"));
        bvhRefractionLabel = new QLabel(bvh);
        bvhRefractionLabel->setObjectName(QStringLiteral("bvhRefractionLabel"));
        bvhRefractionLabel->setFont(font2);

        horizontalLayout_8->addWidget(bvhRefractionLabel);

        bvhRefraction = new QLineEdit(bvh);
        bvhRefraction->setObjectName(QStringLiteral("bvhRefraction"));
        bvhRefraction->setMaximumSize(QSize(100, 20));

        horizontalLayout_8->addWidget(bvhRefraction);


        verticalLayout_5->addLayout(horizontalLayout_8);

        label_9 = new QLabel(bvh);
        label_9->setObjectName(QStringLiteral("label_9"));
        label_9->setFont(font3);
        label_9->setWordWrap(true);

        verticalLayout_5->addWidget(label_9);

        verticalSpacer_3 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_5->addItem(verticalSpacer_3);

        renderModesSettings->addWidget(bvh);
        pssmlt = new QWidget();
        pssmlt->setObjectName(QStringLiteral("pssmlt"));
        verticalLayout_6 = new QVBoxLayout(pssmlt);
        verticalLayout_6->setSpacing(6);
        verticalLayout_6->setContentsMargins(11, 11, 11, 11);
        verticalLayout_6->setObjectName(QStringLiteral("verticalLayout_6"));
        horizontalLayout_11 = new QHBoxLayout();
        horizontalLayout_11->setSpacing(6);
        horizontalLayout_11->setObjectName(QStringLiteral("horizontalLayout_11"));
        pssmltSizeLabel = new QLabel(pssmlt);
        pssmltSizeLabel->setObjectName(QStringLiteral("pssmltSizeLabel"));
        pssmltSizeLabel->setFont(font2);

        horizontalLayout_11->addWidget(pssmltSizeLabel);

        pssmltSize = new QLineEdit(pssmlt);
        pssmltSize->setObjectName(QStringLiteral("pssmltSize"));
        pssmltSize->setMaximumSize(QSize(100, 20));

        horizontalLayout_11->addWidget(pssmltSize);


        verticalLayout_6->addLayout(horizontalLayout_11);

        label_10 = new QLabel(pssmlt);
        label_10->setObjectName(QStringLiteral("label_10"));
        label_10->setWordWrap(true);

        verticalLayout_6->addWidget(label_10);

        line_7 = new QFrame(pssmlt);
        line_7->setObjectName(QStringLiteral("line_7"));
        line_7->setFrameShape(QFrame::HLine);
        line_7->setFrameShadow(QFrame::Sunken);

        verticalLayout_6->addWidget(line_7);

        horizontalLayout_12 = new QHBoxLayout();
        horizontalLayout_12->setSpacing(6);
        horizontalLayout_12->setObjectName(QStringLiteral("horizontalLayout_12"));
        pssmltSamplesLabel = new QLabel(pssmlt);
        pssmltSamplesLabel->setObjectName(QStringLiteral("pssmltSamplesLabel"));
        pssmltSamplesLabel->setFont(font2);

        horizontalLayout_12->addWidget(pssmltSamplesLabel);

        pssmltSamples = new QLineEdit(pssmlt);
        pssmltSamples->setObjectName(QStringLiteral("pssmltSamples"));
        pssmltSamples->setMaximumSize(QSize(100, 20));

        horizontalLayout_12->addWidget(pssmltSamples);


        verticalLayout_6->addLayout(horizontalLayout_12);

        label_11 = new QLabel(pssmlt);
        label_11->setObjectName(QStringLiteral("label_11"));
        label_11->setWordWrap(true);

        verticalLayout_6->addWidget(label_11);

        line_8 = new QFrame(pssmlt);
        line_8->setObjectName(QStringLiteral("line_8"));
        line_8->setFrameShape(QFrame::HLine);
        line_8->setFrameShadow(QFrame::Sunken);

        verticalLayout_6->addWidget(line_8);

        horizontalLayout_13 = new QHBoxLayout();
        horizontalLayout_13->setSpacing(6);
        horizontalLayout_13->setObjectName(QStringLiteral("horizontalLayout_13"));
        pssmltRefractionLabel = new QLabel(pssmlt);
        pssmltRefractionLabel->setObjectName(QStringLiteral("pssmltRefractionLabel"));
        pssmltRefractionLabel->setFont(font2);

        horizontalLayout_13->addWidget(pssmltRefractionLabel);

        pssmltRefraction = new QLineEdit(pssmlt);
        pssmltRefraction->setObjectName(QStringLiteral("pssmltRefraction"));
        pssmltRefraction->setMaximumSize(QSize(100, 20));

        horizontalLayout_13->addWidget(pssmltRefraction);


        verticalLayout_6->addLayout(horizontalLayout_13);

        label_12 = new QLabel(pssmlt);
        label_12->setObjectName(QStringLiteral("label_12"));
        label_12->setWordWrap(true);

        verticalLayout_6->addWidget(label_12);

        verticalSpacer_4 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_6->addItem(verticalSpacer_4);

        renderModesSettings->addWidget(pssmlt);
        smallmedia = new QWidget();
        smallmedia->setObjectName(QStringLiteral("smallmedia"));
        verticalLayout_7 = new QVBoxLayout(smallmedia);
        verticalLayout_7->setSpacing(6);
        verticalLayout_7->setContentsMargins(11, 11, 11, 11);
        verticalLayout_7->setObjectName(QStringLiteral("verticalLayout_7"));
        horizontalLayout_14 = new QHBoxLayout();
        horizontalLayout_14->setSpacing(6);
        horizontalLayout_14->setObjectName(QStringLiteral("horizontalLayout_14"));
        smallmediaSizeLabel = new QLabel(smallmedia);
        smallmediaSizeLabel->setObjectName(QStringLiteral("smallmediaSizeLabel"));
        smallmediaSizeLabel->setFont(font2);

        horizontalLayout_14->addWidget(smallmediaSizeLabel);

        smallmediaSize = new QLineEdit(smallmedia);
        smallmediaSize->setObjectName(QStringLiteral("smallmediaSize"));
        smallmediaSize->setMaximumSize(QSize(100, 20));

        horizontalLayout_14->addWidget(smallmediaSize);


        verticalLayout_7->addLayout(horizontalLayout_14);

        label_13 = new QLabel(smallmedia);
        label_13->setObjectName(QStringLiteral("label_13"));
        label_13->setWordWrap(true);

        verticalLayout_7->addWidget(label_13);

        line_9 = new QFrame(smallmedia);
        line_9->setObjectName(QStringLiteral("line_9"));
        line_9->setFrameShape(QFrame::HLine);
        line_9->setFrameShadow(QFrame::Sunken);

        verticalLayout_7->addWidget(line_9);

        horizontalLayout_15 = new QHBoxLayout();
        horizontalLayout_15->setSpacing(6);
        horizontalLayout_15->setObjectName(QStringLiteral("horizontalLayout_15"));
        smallmediaSamplesLabel = new QLabel(smallmedia);
        smallmediaSamplesLabel->setObjectName(QStringLiteral("smallmediaSamplesLabel"));
        smallmediaSamplesLabel->setFont(font2);

        horizontalLayout_15->addWidget(smallmediaSamplesLabel);

        smallmediaSamples = new QLineEdit(smallmedia);
        smallmediaSamples->setObjectName(QStringLiteral("smallmediaSamples"));
        smallmediaSamples->setMaximumSize(QSize(100, 20));

        horizontalLayout_15->addWidget(smallmediaSamples);


        verticalLayout_7->addLayout(horizontalLayout_15);

        label_14 = new QLabel(smallmedia);
        label_14->setObjectName(QStringLiteral("label_14"));
        label_14->setWordWrap(true);

        verticalLayout_7->addWidget(label_14);

        line_10 = new QFrame(smallmedia);
        line_10->setObjectName(QStringLiteral("line_10"));
        line_10->setFrameShape(QFrame::HLine);
        line_10->setFrameShadow(QFrame::Sunken);

        verticalLayout_7->addWidget(line_10);

        horizontalLayout_16 = new QHBoxLayout();
        horizontalLayout_16->setSpacing(6);
        horizontalLayout_16->setObjectName(QStringLiteral("horizontalLayout_16"));
        smallmediaRefractionLabel = new QLabel(smallmedia);
        smallmediaRefractionLabel->setObjectName(QStringLiteral("smallmediaRefractionLabel"));
        smallmediaRefractionLabel->setFont(font2);

        horizontalLayout_16->addWidget(smallmediaRefractionLabel);

        smallmediaRefraction = new QLineEdit(smallmedia);
        smallmediaRefraction->setObjectName(QStringLiteral("smallmediaRefraction"));
        smallmediaRefraction->setMaximumSize(QSize(100, 20));

        horizontalLayout_16->addWidget(smallmediaRefraction);


        verticalLayout_7->addLayout(horizontalLayout_16);

        label_15 = new QLabel(smallmedia);
        label_15->setObjectName(QStringLiteral("label_15"));
        label_15->setWordWrap(true);

        verticalLayout_7->addWidget(label_15);

        verticalSpacer_5 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_7->addItem(verticalSpacer_5);

        renderModesSettings->addWidget(smallmedia);
        ppm = new QWidget();
        ppm->setObjectName(QStringLiteral("ppm"));
        verticalLayout_8 = new QVBoxLayout(ppm);
        verticalLayout_8->setSpacing(6);
        verticalLayout_8->setContentsMargins(11, 11, 11, 11);
        verticalLayout_8->setObjectName(QStringLiteral("verticalLayout_8"));
        horizontalLayout_17 = new QHBoxLayout();
        horizontalLayout_17->setSpacing(6);
        horizontalLayout_17->setObjectName(QStringLiteral("horizontalLayout_17"));
        ppmSizeLabel = new QLabel(ppm);
        ppmSizeLabel->setObjectName(QStringLiteral("ppmSizeLabel"));
        ppmSizeLabel->setFont(font2);

        horizontalLayout_17->addWidget(ppmSizeLabel);

        ppmSize = new QLineEdit(ppm);
        ppmSize->setObjectName(QStringLiteral("ppmSize"));
        ppmSize->setMaximumSize(QSize(100, 20));

        horizontalLayout_17->addWidget(ppmSize);


        verticalLayout_8->addLayout(horizontalLayout_17);

        label_16 = new QLabel(ppm);
        label_16->setObjectName(QStringLiteral("label_16"));
        label_16->setWordWrap(true);

        verticalLayout_8->addWidget(label_16);

        line_11 = new QFrame(ppm);
        line_11->setObjectName(QStringLiteral("line_11"));
        line_11->setFrameShape(QFrame::HLine);
        line_11->setFrameShadow(QFrame::Sunken);

        verticalLayout_8->addWidget(line_11);

        horizontalLayout_18 = new QHBoxLayout();
        horizontalLayout_18->setSpacing(6);
        horizontalLayout_18->setObjectName(QStringLiteral("horizontalLayout_18"));
        ppmSamplesLabel = new QLabel(ppm);
        ppmSamplesLabel->setObjectName(QStringLiteral("ppmSamplesLabel"));
        ppmSamplesLabel->setFont(font2);

        horizontalLayout_18->addWidget(ppmSamplesLabel);

        ppmSamples = new QLineEdit(ppm);
        ppmSamples->setObjectName(QStringLiteral("ppmSamples"));
        ppmSamples->setMaximumSize(QSize(100, 20));

        horizontalLayout_18->addWidget(ppmSamples);


        verticalLayout_8->addLayout(horizontalLayout_18);

        label_17 = new QLabel(ppm);
        label_17->setObjectName(QStringLiteral("label_17"));
        label_17->setWordWrap(true);

        verticalLayout_8->addWidget(label_17);

        line_12 = new QFrame(ppm);
        line_12->setObjectName(QStringLiteral("line_12"));
        line_12->setFrameShape(QFrame::HLine);
        line_12->setFrameShadow(QFrame::Sunken);

        verticalLayout_8->addWidget(line_12);

        horizontalLayout_19 = new QHBoxLayout();
        horizontalLayout_19->setSpacing(6);
        horizontalLayout_19->setObjectName(QStringLiteral("horizontalLayout_19"));
        ppmRefractionLabel = new QLabel(ppm);
        ppmRefractionLabel->setObjectName(QStringLiteral("ppmRefractionLabel"));
        ppmRefractionLabel->setFont(font2);

        horizontalLayout_19->addWidget(ppmRefractionLabel);

        ppmRefraction = new QLineEdit(ppm);
        ppmRefraction->setObjectName(QStringLiteral("ppmRefraction"));
        ppmRefraction->setMaximumSize(QSize(100, 20));

        horizontalLayout_19->addWidget(ppmRefraction);


        verticalLayout_8->addLayout(horizontalLayout_19);

        label_18 = new QLabel(ppm);
        label_18->setObjectName(QStringLiteral("label_18"));
        label_18->setWordWrap(true);

        verticalLayout_8->addWidget(label_18);

        verticalSpacer_6 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_8->addItem(verticalSpacer_6);

        renderModesSettings->addWidget(ppm);
        vrl = new QWidget();
        vrl->setObjectName(QStringLiteral("vrl"));
        verticalLayout_2 = new QVBoxLayout(vrl);
        verticalLayout_2->setSpacing(0);
        verticalLayout_2->setContentsMargins(11, 11, 11, 11);
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        verticalLayout_2->setContentsMargins(0, 0, 0, 0);
        scrollArea = new QScrollArea(vrl);
        scrollArea->setObjectName(QStringLiteral("scrollArea"));
        scrollArea->setFrameShape(QFrame::NoFrame);
        scrollArea->setFrameShadow(QFrame::Plain);
        scrollArea->setLineWidth(0);
        scrollArea->setWidgetResizable(true);
        scrollAreaWidgetContents = new QWidget();
        scrollAreaWidgetContents->setObjectName(QStringLiteral("scrollAreaWidgetContents"));
        scrollAreaWidgetContents->setGeometry(QRect(0, 0, 188, 604));
        scrollAreaWidgetContents->setFont(font1);
        scrollAreaWidgetContents->setAutoFillBackground(true);
        verticalLayout_9 = new QVBoxLayout(scrollAreaWidgetContents);
        verticalLayout_9->setSpacing(6);
        verticalLayout_9->setContentsMargins(11, 11, 11, 11);
        verticalLayout_9->setObjectName(QStringLiteral("verticalLayout_9"));
        verticalLayout_9->setContentsMargins(9, 9, 9, 9);
        horizontalLayout_22 = new QHBoxLayout();
        horizontalLayout_22->setSpacing(6);
        horizontalLayout_22->setObjectName(QStringLiteral("horizontalLayout_22"));
        vrlSizeLabel = new QLabel(scrollAreaWidgetContents);
        vrlSizeLabel->setObjectName(QStringLiteral("vrlSizeLabel"));
        vrlSizeLabel->setFont(font2);

        horizontalLayout_22->addWidget(vrlSizeLabel);

        vrlSize = new QLineEdit(scrollAreaWidgetContents);
        vrlSize->setObjectName(QStringLiteral("vrlSize"));
        vrlSize->setMaximumSize(QSize(100, 20));

        horizontalLayout_22->addWidget(vrlSize);


        verticalLayout_9->addLayout(horizontalLayout_22);

        label_19 = new QLabel(scrollAreaWidgetContents);
        label_19->setObjectName(QStringLiteral("label_19"));
        label_19->setWordWrap(true);

        verticalLayout_9->addWidget(label_19);

        line_13 = new QFrame(scrollAreaWidgetContents);
        line_13->setObjectName(QStringLiteral("line_13"));
        line_13->setFrameShape(QFrame::HLine);
        line_13->setFrameShadow(QFrame::Sunken);

        verticalLayout_9->addWidget(line_13);

        horizontalLayout_20 = new QHBoxLayout();
        horizontalLayout_20->setSpacing(6);
        horizontalLayout_20->setObjectName(QStringLiteral("horizontalLayout_20"));
        vrlSamplesLabel = new QLabel(scrollAreaWidgetContents);
        vrlSamplesLabel->setObjectName(QStringLiteral("vrlSamplesLabel"));
        vrlSamplesLabel->setFont(font2);

        horizontalLayout_20->addWidget(vrlSamplesLabel);

        vrlSamples = new QLineEdit(scrollAreaWidgetContents);
        vrlSamples->setObjectName(QStringLiteral("vrlSamples"));
        vrlSamples->setMaximumSize(QSize(100, 20));

        horizontalLayout_20->addWidget(vrlSamples);


        verticalLayout_9->addLayout(horizontalLayout_20);

        label_20 = new QLabel(scrollAreaWidgetContents);
        label_20->setObjectName(QStringLiteral("label_20"));
        label_20->setWordWrap(true);

        verticalLayout_9->addWidget(label_20);

        line_14 = new QFrame(scrollAreaWidgetContents);
        line_14->setObjectName(QStringLiteral("line_14"));
        line_14->setFrameShape(QFrame::HLine);
        line_14->setFrameShadow(QFrame::Sunken);

        verticalLayout_9->addWidget(line_14);

        horizontalLayout_23 = new QHBoxLayout();
        horizontalLayout_23->setSpacing(6);
        horizontalLayout_23->setObjectName(QStringLiteral("horizontalLayout_23"));
        vrlVRLPSLabel = new QLabel(scrollAreaWidgetContents);
        vrlVRLPSLabel->setObjectName(QStringLiteral("vrlVRLPSLabel"));
        vrlVRLPSLabel->setFont(font2);

        horizontalLayout_23->addWidget(vrlVRLPSLabel);

        vrlVRLPS = new QLineEdit(scrollAreaWidgetContents);
        vrlVRLPS->setObjectName(QStringLiteral("vrlVRLPS"));
        vrlVRLPS->setMaximumSize(QSize(100, 20));

        horizontalLayout_23->addWidget(vrlVRLPS);


        verticalLayout_9->addLayout(horizontalLayout_23);

        label_23 = new QLabel(scrollAreaWidgetContents);
        label_23->setObjectName(QStringLiteral("label_23"));
        label_23->setWordWrap(true);

        verticalLayout_9->addWidget(label_23);

        line_16 = new QFrame(scrollAreaWidgetContents);
        line_16->setObjectName(QStringLiteral("line_16"));
        line_16->setFrameShape(QFrame::HLine);
        line_16->setFrameShadow(QFrame::Sunken);

        verticalLayout_9->addWidget(line_16);

        horizontalLayout_26 = new QHBoxLayout();
        horizontalLayout_26->setSpacing(6);
        horizontalLayout_26->setObjectName(QStringLiteral("horizontalLayout_26"));
        vrlBouncesLabel = new QLabel(scrollAreaWidgetContents);
        vrlBouncesLabel->setObjectName(QStringLiteral("vrlBouncesLabel"));
        vrlBouncesLabel->setFont(font2);

        horizontalLayout_26->addWidget(vrlBouncesLabel);

        vrlBounces = new QLineEdit(scrollAreaWidgetContents);
        vrlBounces->setObjectName(QStringLiteral("vrlBounces"));
        vrlBounces->setMaximumSize(QSize(100, 16777215));

        horizontalLayout_26->addWidget(vrlBounces);


        verticalLayout_9->addLayout(horizontalLayout_26);

        label_27 = new QLabel(scrollAreaWidgetContents);
        label_27->setObjectName(QStringLiteral("label_27"));
        label_27->setWordWrap(true);

        verticalLayout_9->addWidget(label_27);

        line_17 = new QFrame(scrollAreaWidgetContents);
        line_17->setObjectName(QStringLiteral("line_17"));
        line_17->setFrameShape(QFrame::HLine);
        line_17->setFrameShadow(QFrame::Sunken);

        verticalLayout_9->addWidget(line_17);

        horizontalLayout_27 = new QHBoxLayout();
        horizontalLayout_27->setSpacing(6);
        horizontalLayout_27->setObjectName(QStringLiteral("horizontalLayout_27"));
        vrlRefractionLabel = new QLabel(scrollAreaWidgetContents);
        vrlRefractionLabel->setObjectName(QStringLiteral("vrlRefractionLabel"));
        vrlRefractionLabel->setFont(font2);

        horizontalLayout_27->addWidget(vrlRefractionLabel);

        vrlRefraction = new QLineEdit(scrollAreaWidgetContents);
        vrlRefraction->setObjectName(QStringLiteral("vrlRefraction"));
        vrlRefraction->setMaximumSize(QSize(100, 16777215));

        horizontalLayout_27->addWidget(vrlRefraction);


        verticalLayout_9->addLayout(horizontalLayout_27);

        label_29 = new QLabel(scrollAreaWidgetContents);
        label_29->setObjectName(QStringLiteral("label_29"));
        label_29->setWordWrap(true);

        verticalLayout_9->addWidget(label_29);

        line_19 = new QFrame(scrollAreaWidgetContents);
        line_19->setObjectName(QStringLiteral("line_19"));
        line_19->setFrameShape(QFrame::HLine);
        line_19->setFrameShadow(QFrame::Sunken);

        verticalLayout_9->addWidget(line_19);

        horizontalLayout_25 = new QHBoxLayout();
        horizontalLayout_25->setSpacing(6);
        horizontalLayout_25->setObjectName(QStringLiteral("horizontalLayout_25"));
        vrlMediumLabel = new QLabel(scrollAreaWidgetContents);
        vrlMediumLabel->setObjectName(QStringLiteral("vrlMediumLabel"));
        vrlMediumLabel->setFont(font2);

        horizontalLayout_25->addWidget(vrlMediumLabel);

        vrlSigmaA = new QLineEdit(scrollAreaWidgetContents);
        vrlSigmaA->setObjectName(QStringLiteral("vrlSigmaA"));
        vrlSigmaA->setMaximumSize(QSize(35, 16777215));

        horizontalLayout_25->addWidget(vrlSigmaA);

        vrlSigmaS = new QLineEdit(scrollAreaWidgetContents);
        vrlSigmaS->setObjectName(QStringLiteral("vrlSigmaS"));
        vrlSigmaS->setMaximumSize(QSize(35, 16777215));

        horizontalLayout_25->addWidget(vrlSigmaS);

        vrlG = new QLineEdit(scrollAreaWidgetContents);
        vrlG->setObjectName(QStringLiteral("vrlG"));
        vrlG->setMaximumSize(QSize(35, 16777215));

        horizontalLayout_25->addWidget(vrlG);


        verticalLayout_9->addLayout(horizontalLayout_25);

        label_25 = new QLabel(scrollAreaWidgetContents);
        label_25->setObjectName(QStringLiteral("label_25"));
        label_25->setWordWrap(true);

        verticalLayout_9->addWidget(label_25);

        line_21 = new QFrame(scrollAreaWidgetContents);
        line_21->setObjectName(QStringLiteral("line_21"));
        line_21->setFrameShape(QFrame::HLine);
        line_21->setFrameShadow(QFrame::Sunken);

        verticalLayout_9->addWidget(line_21);

        horizontalLayout_28 = new QHBoxLayout();
        horizontalLayout_28->setSpacing(6);
        horizontalLayout_28->setObjectName(QStringLiteral("horizontalLayout_28"));
        vrlIntensityLabel = new QLabel(scrollAreaWidgetContents);
        vrlIntensityLabel->setObjectName(QStringLiteral("vrlIntensityLabel"));
        vrlIntensityLabel->setFont(font2);

        horizontalLayout_28->addWidget(vrlIntensityLabel);

        vrlIntensity = new QLineEdit(scrollAreaWidgetContents);
        vrlIntensity->setObjectName(QStringLiteral("vrlIntensity"));
        vrlIntensity->setMaximumSize(QSize(100, 16777215));

        horizontalLayout_28->addWidget(vrlIntensity);


        verticalLayout_9->addLayout(horizontalLayout_28);

        label_26 = new QLabel(scrollAreaWidgetContents);
        label_26->setObjectName(QStringLiteral("label_26"));
        label_26->setWordWrap(true);

        verticalLayout_9->addWidget(label_26);

        verticalSpacer_7 = new QSpacerItem(20, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_9->addItem(verticalSpacer_7);

        line_20 = new QFrame(scrollAreaWidgetContents);
        line_20->setObjectName(QStringLiteral("line_20"));
        line_20->setFrameShape(QFrame::HLine);
        line_20->setFrameShadow(QFrame::Sunken);

        verticalLayout_9->addWidget(line_20);

        vrlMediumRadiance = new QCheckBox(scrollAreaWidgetContents);
        vrlMediumRadiance->setObjectName(QStringLiteral("vrlMediumRadiance"));
        vrlMediumRadiance->setFont(font2);
        vrlMediumRadiance->setChecked(true);

        verticalLayout_9->addWidget(vrlMediumRadiance);

        vrlSurfaceRadiance = new QCheckBox(scrollAreaWidgetContents);
        vrlSurfaceRadiance->setObjectName(QStringLiteral("vrlSurfaceRadiance"));
        vrlSurfaceRadiance->setFont(font2);
        vrlSurfaceRadiance->setChecked(true);

        verticalLayout_9->addWidget(vrlSurfaceRadiance);

        line_15 = new QFrame(scrollAreaWidgetContents);
        line_15->setObjectName(QStringLiteral("line_15"));
        line_15->setFrameShape(QFrame::HLine);
        line_15->setFrameShadow(QFrame::Sunken);

        verticalLayout_9->addWidget(line_15);

        horizontalLayout_21 = new QHBoxLayout();
        horizontalLayout_21->setSpacing(6);
        horizontalLayout_21->setObjectName(QStringLiteral("horizontalLayout_21"));
        label_21 = new QLabel(scrollAreaWidgetContents);
        label_21->setObjectName(QStringLiteral("label_21"));
        label_21->setFont(font2);

        horizontalLayout_21->addWidget(label_21);

        vrlSampling = new QComboBox(scrollAreaWidgetContents);
        vrlSampling->setObjectName(QStringLiteral("vrlSampling"));
        vrlSampling->setMinimumSize(QSize(0, 0));
        vrlSampling->setMaximumSize(QSize(100, 16777215));

        horizontalLayout_21->addWidget(vrlSampling);


        verticalLayout_9->addLayout(horizontalLayout_21);

        label_22 = new QLabel(scrollAreaWidgetContents);
        label_22->setObjectName(QStringLiteral("label_22"));
        label_22->setWordWrap(true);

        verticalLayout_9->addWidget(label_22);

        line_18 = new QFrame(scrollAreaWidgetContents);
        line_18->setObjectName(QStringLiteral("line_18"));
        line_18->setFrameShape(QFrame::HLine);
        line_18->setFrameShadow(QFrame::Sunken);

        verticalLayout_9->addWidget(line_18);

        horizontalLayout_24 = new QHBoxLayout();
        horizontalLayout_24->setSpacing(6);
        horizontalLayout_24->setObjectName(QStringLiteral("horizontalLayout_24"));
        label_24 = new QLabel(scrollAreaWidgetContents);
        label_24->setObjectName(QStringLiteral("label_24"));
        label_24->setFont(font2);

        horizontalLayout_24->addWidget(label_24);

        vrlScene = new QComboBox(scrollAreaWidgetContents);
        vrlScene->setObjectName(QStringLiteral("vrlScene"));
        vrlScene->setMaximumSize(QSize(100, 16777215));

        horizontalLayout_24->addWidget(vrlScene);


        verticalLayout_9->addLayout(horizontalLayout_24);

        scrollArea->setWidget(scrollAreaWidgetContents);

        verticalLayout_2->addWidget(scrollArea);

        renderModesSettings->addWidget(vrl);
        tests = new QWidget();
        tests->setObjectName(QStringLiteral("tests"));
        formLayout_5 = new QFormLayout(tests);
        formLayout_5->setSpacing(6);
        formLayout_5->setContentsMargins(11, 11, 11, 11);
        formLayout_5->setObjectName(QStringLiteral("formLayout_5"));
        fixedTest = new QRadioButton(tests);
        fixedTest->setObjectName(QStringLiteral("fixedTest"));
        fixedTest->setMaximumSize(QSize(300, 16777215));
        fixedTest->setChecked(true);

        formLayout_5->setWidget(0, QFormLayout::LabelRole, fixedTest);

        painterlyTest = new QRadioButton(tests);
        painterlyTest->setObjectName(QStringLiteral("painterlyTest"));
        painterlyTest->setMaximumSize(QSize(300, 16777215));
        painterlyTest->setChecked(false);

        formLayout_5->setWidget(1, QFormLayout::LabelRole, painterlyTest);

        bvhTest = new QRadioButton(tests);
        bvhTest->setObjectName(QStringLiteral("bvhTest"));
        bvhTest->setMaximumSize(QSize(300, 16777215));
        bvhTest->setChecked(false);

        formLayout_5->setWidget(2, QFormLayout::LabelRole, bvhTest);

        pssmltTest = new QRadioButton(tests);
        pssmltTest->setObjectName(QStringLiteral("pssmltTest"));
        pssmltTest->setMaximumSize(QSize(300, 16777215));
        pssmltTest->setChecked(false);

        formLayout_5->setWidget(3, QFormLayout::LabelRole, pssmltTest);

        smallmediaTest = new QRadioButton(tests);
        smallmediaTest->setObjectName(QStringLiteral("smallmediaTest"));
        smallmediaTest->setMaximumSize(QSize(300, 16777215));
        smallmediaTest->setChecked(false);

        formLayout_5->setWidget(4, QFormLayout::LabelRole, smallmediaTest);

        ppmTest = new QRadioButton(tests);
        ppmTest->setObjectName(QStringLiteral("ppmTest"));
        ppmTest->setMaximumSize(QSize(300, 16777215));
        ppmTest->setChecked(false);

        formLayout_5->setWidget(5, QFormLayout::LabelRole, ppmTest);

        vrlTest = new QRadioButton(tests);
        vrlTest->setObjectName(QStringLiteral("vrlTest"));

        formLayout_5->setWidget(6, QFormLayout::LabelRole, vrlTest);

        renderModesSettings->addWidget(tests);

        modes->addWidget(renderModesSettings);

        timeInfo = new QTextBrowser(frame);
        timeInfo->setObjectName(QStringLiteral("timeInfo"));
        sizePolicy.setHeightForWidth(timeInfo->sizePolicy().hasHeightForWidth());
        timeInfo->setSizePolicy(sizePolicy);
        timeInfo->setMinimumSize(QSize(0, 0));
        timeInfo->setMaximumSize(QSize(300, 100));
        timeInfo->setAutoFillBackground(true);
        timeInfo->setFrameShape(QFrame::NoFrame);
        timeInfo->setFrameShadow(QFrame::Sunken);
        timeInfo->setLineWidth(1);

        modes->addWidget(timeInfo);

        progressBar = new QProgressBar(frame);
        progressBar->setObjectName(QStringLiteral("progressBar"));
        QSizePolicy sizePolicy3(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(progressBar->sizePolicy().hasHeightForWidth());
        progressBar->setSizePolicy(sizePolicy3);
        progressBar->setMinimumSize(QSize(150, 28));
        progressBar->setMaximumSize(QSize(300, 16777215));
        progressBar->setValue(0);
        progressBar->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignVCenter);
        progressBar->setTextVisible(false);
        progressBar->setInvertedAppearance(false);

        modes->addWidget(progressBar);


        horizontalLayout->addWidget(frame);

        image = new Image(centralWidget);
        image->setObjectName(QStringLiteral("image"));
        sizePolicy2.setHeightForWidth(image->sizePolicy().hasHeightForWidth());
        image->setSizePolicy(sizePolicy2);
        image->setMinimumSize(QSize(512, 512));

        horizontalLayout->addWidget(image);


        verticalLayout->addLayout(horizontalLayout);

        horizontalLayoutBottom = new QHBoxLayout();
        horizontalLayoutBottom->setSpacing(6);
        horizontalLayoutBottom->setObjectName(QStringLiteral("horizontalLayoutBottom"));
        renderInfo = new QTextBrowser(centralWidget);
        renderInfo->setObjectName(QStringLiteral("renderInfo"));
        sizePolicy3.setHeightForWidth(renderInfo->sizePolicy().hasHeightForWidth());
        renderInfo->setSizePolicy(sizePolicy3);
        renderInfo->setMaximumSize(QSize(16777215, 30));
        renderInfo->setAutoFillBackground(false);
        renderInfo->setFrameShape(QFrame::NoFrame);
        renderInfo->setFrameShadow(QFrame::Sunken);
        renderInfo->setLineWidth(1);
        renderInfo->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
        renderInfo->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);

        horizontalLayoutBottom->addWidget(renderInfo);

        renderButton = new QPushButton(centralWidget);
        renderButton->setObjectName(QStringLiteral("renderButton"));
        QSizePolicy sizePolicy4(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy4.setHorizontalStretch(0);
        sizePolicy4.setVerticalStretch(0);
        sizePolicy4.setHeightForWidth(renderButton->sizePolicy().hasHeightForWidth());
        renderButton->setSizePolicy(sizePolicy4);
        renderButton->setMinimumSize(QSize(120, 30));
        renderButton->setToolTipDuration(-1);
        renderButton->setIconSize(QSize(16, 16));
        renderButton->setAutoRepeatInterval(100);
        renderButton->setAutoDefault(false);

        horizontalLayoutBottom->addWidget(renderButton);


        verticalLayout->addLayout(horizontalLayoutBottom);

        MainWindow->setCentralWidget(centralWidget);

        retranslateUi(MainWindow);

        renderModesComboBox->setCurrentIndex(0);
        renderModesSettings->setCurrentIndex(0);
        vrlSampling->setCurrentIndex(4);
        renderButton->setDefault(false);


        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "Smallpaint", 0));
        renderModesComboBox->clear();
        renderModesComboBox->insertItems(0, QStringList()
         << QApplication::translate("MainWindow", "Path Tracing (Fixed)", 0)
         << QApplication::translate("MainWindow", "Path Tracing (Painterly)", 0)
         << QApplication::translate("MainWindow", "Bounding Volume Hierarchy", 0)
         << QApplication::translate("MainWindow", "Primary Sample Space MLT", 0)
         << QApplication::translate("MainWindow", "Volumetric Path Tracing", 0)
         << QApplication::translate("MainWindow", "Progressive Photon Mapping", 0)
         << QApplication::translate("MainWindow", "Virtual Ray Lights", 0)
         << QApplication::translate("MainWindow", "Tests", 0)
        );
#ifndef QT_NO_TOOLTIP
        fixedSizeLabel->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        fixedSizeLabel->setText(QApplication::translate("MainWindow", "Width / Height", 0));
#ifndef QT_NO_TOOLTIP
        fixedSize->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        fixedSize->setText(QApplication::translate("MainWindow", "512", 0));
        label->setText(QApplication::translate("MainWindow", "The width and height of the image. Entering 400 creates a 400x400 image.", 0));
        fixedSamplesLabel->setText(QApplication::translate("MainWindow", "Samples", 0));
#ifndef QT_NO_TOOLTIP
        fixedSamples->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        fixedSamples->setText(QApplication::translate("MainWindow", "50", 0));
        label_2->setText(QApplication::translate("MainWindow", "The number of samples per pixel that are used to render the image.", 0));
        fixedRefractionLabel->setText(QApplication::translate("MainWindow", "Refraction Index", 0));
#ifndef QT_NO_TOOLTIP
        fixedRefraction->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        fixedRefraction->setText(QApplication::translate("MainWindow", "1.5", 0));
        label_3->setText(QApplication::translate("MainWindow", "The refraction index for refractive spheres. (Water 1.3; Diamond 2.4)", 0));
        painterlySizeLabel->setText(QApplication::translate("MainWindow", "Width / Height", 0));
#ifndef QT_NO_TOOLTIP
        painterlySize->setToolTip(QApplication::translate("MainWindow", "<html><head/><body><p><br/></p></body></html>", 0));
#endif // QT_NO_TOOLTIP
        label_4->setText(QApplication::translate("MainWindow", "The width and height of the image. Entering 400 creates a 400x400 image.", 0));
        painterlySamplesLabel->setText(QApplication::translate("MainWindow", "Samples", 0));
#ifndef QT_NO_TOOLTIP
        painterlySamples->setToolTip(QApplication::translate("MainWindow", "<html><head/><body><p><br/></p></body></html>", 0));
#endif // QT_NO_TOOLTIP
        label_5->setText(QApplication::translate("MainWindow", "The number of samples per pixel that are used to render the image.", 0));
        painterlyRefractionLabel->setText(QApplication::translate("MainWindow", "Refraction Index", 0));
#ifndef QT_NO_TOOLTIP
        painterlyRefraction->setToolTip(QApplication::translate("MainWindow", "<html><head/><body><p><br/></p></body></html>", 0));
#endif // QT_NO_TOOLTIP
        label_6->setText(QApplication::translate("MainWindow", "The refraction index for refractive spheres. (Water 1.3; Diamond 2.4)", 0));
        bvhSizeLabel->setText(QApplication::translate("MainWindow", "Width / Height", 0));
#ifndef QT_NO_TOOLTIP
        bvhSize->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        label_7->setText(QApplication::translate("MainWindow", "The width and height of the image. Entering 400 creates a 400x400 image.", 0));
        bvhSamplesLabel->setText(QApplication::translate("MainWindow", "Samples", 0));
#ifndef QT_NO_TOOLTIP
        bvhSamples->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        label_8->setText(QApplication::translate("MainWindow", "The number of samples per pixel that are used to render the image.", 0));
        bvhRefractionLabel->setText(QApplication::translate("MainWindow", "Refraction Index", 0));
#ifndef QT_NO_TOOLTIP
        bvhRefraction->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        label_9->setText(QApplication::translate("MainWindow", "The refraction index for refractive spheres. (Water 1.3; Diamond 2.4)", 0));
        pssmltSizeLabel->setText(QApplication::translate("MainWindow", "Width / Height", 0));
#ifndef QT_NO_TOOLTIP
        pssmltSize->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        label_10->setText(QApplication::translate("MainWindow", "The width and height of the image. Entering 400 creates a 400x400 image.", 0));
        pssmltSamplesLabel->setText(QApplication::translate("MainWindow", "Samples", 0));
#ifndef QT_NO_TOOLTIP
        pssmltSamples->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        label_11->setText(QApplication::translate("MainWindow", "The number of samples per pixel that are used to render the image.", 0));
        pssmltRefractionLabel->setText(QApplication::translate("MainWindow", "Refraction Index", 0));
#ifndef QT_NO_TOOLTIP
        pssmltRefraction->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        label_12->setText(QApplication::translate("MainWindow", "The refraction index for refractive spheres. (Water 1.3; Diamond 2.4)", 0));
        smallmediaSizeLabel->setText(QApplication::translate("MainWindow", "Width / Height", 0));
#ifndef QT_NO_TOOLTIP
        smallmediaSize->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        label_13->setText(QApplication::translate("MainWindow", "The width and height of the image. Entering 400 creates a 400x400 image.", 0));
        smallmediaSamplesLabel->setText(QApplication::translate("MainWindow", "Samples", 0));
#ifndef QT_NO_TOOLTIP
        smallmediaSamples->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        label_14->setText(QApplication::translate("MainWindow", "The number of samples per pixel that are used to render the image.", 0));
        smallmediaRefractionLabel->setText(QApplication::translate("MainWindow", "Refraction Index", 0));
#ifndef QT_NO_TOOLTIP
        smallmediaRefraction->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        label_15->setText(QApplication::translate("MainWindow", "The refraction index for refractive spheres. (Water 1.3; Diamond 2.4)", 0));
        ppmSizeLabel->setText(QApplication::translate("MainWindow", "Width / Height", 0));
#ifndef QT_NO_TOOLTIP
        ppmSize->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        label_16->setText(QApplication::translate("MainWindow", "The width and height of the image. Entering 400 creates a 400x400 image.", 0));
        ppmSamplesLabel->setText(QApplication::translate("MainWindow", "Photons (*1k)", 0));
#ifndef QT_NO_TOOLTIP
        ppmSamples->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        label_17->setText(QApplication::translate("MainWindow", "The number of photons that are used to render the image. The entered number will be multiplied by 1000.", 0));
        ppmRefractionLabel->setText(QApplication::translate("MainWindow", "Refraction Index", 0));
#ifndef QT_NO_TOOLTIP
        ppmRefraction->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        label_18->setText(QApplication::translate("MainWindow", "The refraction index for refractive spheres. (Water 1.3; Diamond 2.4)", 0));
        vrlSizeLabel->setText(QApplication::translate("MainWindow", "Width / Height", 0));
#ifndef QT_NO_TOOLTIP
        vrlSize->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        label_19->setText(QApplication::translate("MainWindow", "The width and height of the image. Entering 400 creates a 400x400 image.", 0));
        vrlSamplesLabel->setText(QApplication::translate("MainWindow", "Samples", 0));
#ifndef QT_NO_TOOLTIP
        vrlSamples->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        label_20->setText(QApplication::translate("MainWindow", "The number of samples per pixel that are used to render the image.", 0));
        vrlVRLPSLabel->setText(QApplication::translate("MainWindow", "VRL per Sample", 0));
        label_23->setText(QApplication::translate("MainWindow", "The number of Virtual Ray Lights per sample that are used to render the image.", 0));
        vrlBouncesLabel->setText(QApplication::translate("MainWindow", "Bounces", 0));
        label_27->setText(QApplication::translate("MainWindow", "The number of allowed bounces for a light ray.", 0));
        vrlRefractionLabel->setText(QApplication::translate("MainWindow", "Refraction Index", 0));
        label_29->setText(QApplication::translate("MainWindow", "The refraction index for refractive spheres. (Water 1.3; Diamond 2.4)", 0));
        vrlMediumLabel->setText(QApplication::translate("MainWindow", "Medium", 0));
        label_25->setText(QApplication::translate("MainWindow", "The medium properties. Absorption coefficient (sigma_a), scattering coefficient (sigma_s) and Henyey-Greensteen g-value.", 0));
        vrlIntensityLabel->setText(QApplication::translate("MainWindow", "Light Intensity", 0));
        label_26->setText(QApplication::translate("MainWindow", "The intensity of the light.", 0));
        vrlMediumRadiance->setText(QApplication::translate("MainWindow", "Medium Radiance", 0));
        vrlSurfaceRadiance->setText(QApplication::translate("MainWindow", "Surface Radiance", 0));
        label_21->setText(QApplication::translate("MainWindow", "Sampling", 0));
        vrlSampling->clear();
        vrlSampling->insertItems(0, QStringList()
         << QApplication::translate("MainWindow", "uniform / uniform", 0)
         << QApplication::translate("MainWindow", "exponential / exponential", 0)
         << QApplication::translate("MainWindow", "uniform / equi-angular", 0)
         << QApplication::translate("MainWindow", "simple joint distribution", 0)
         << QApplication::translate("MainWindow", "advanced joint distribution", 0)
        );
        label_22->setText(QApplication::translate("MainWindow", "The importance sampling strategy for the light rays / camera rays.", 0));
        label_24->setText(QApplication::translate("MainWindow", "Scene", 0));
        vrlScene->clear();
        vrlScene->insertItems(0, QStringList()
         << QApplication::translate("MainWindow", "Cornell Box", 0)
         << QApplication::translate("MainWindow", "Bidirectional VRL", 0)
         << QApplication::translate("MainWindow", "Unidirectional VRL", 0)
         << QApplication::translate("MainWindow", "Medium Radiance Test", 0)
         << QApplication::translate("MainWindow", "Surface Radiance Test", 0)
         << QApplication::translate("MainWindow", "Medium Cube", 0)
        );
        fixedTest->setText(QApplication::translate("MainWindow", "Path Tracing (Fixed)", 0));
        painterlyTest->setText(QApplication::translate("MainWindow", "Path Tracing (Painterly)", 0));
        bvhTest->setText(QApplication::translate("MainWindow", "Bounding Volume Hierarchy", 0));
        pssmltTest->setText(QApplication::translate("MainWindow", "Primary Sample Space MLT", 0));
        smallmediaTest->setText(QApplication::translate("MainWindow", "Volumetric Path Tracing", 0));
        ppmTest->setText(QApplication::translate("MainWindow", "Progressive Photon Mapping", 0));
        vrlTest->setText(QApplication::translate("MainWindow", "Virtual Ray Lights", 0));
        renderButton->setText(QApplication::translate("MainWindow", "Render", 0));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
