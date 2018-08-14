/********************************************************************************
** Form generated from reading UI file 'testwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.7.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_TESTWINDOW_H
#define UI_TESTWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_TestWindow
{
public:
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout;
    QLabel *img1;
    QLabel *img2;

    void setupUi(QDialog *TestWindow)
    {
        if (TestWindow->objectName().isEmpty())
            TestWindow->setObjectName(QStringLiteral("TestWindow"));
        TestWindow->resize(1024, 512);
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(TestWindow->sizePolicy().hasHeightForWidth());
        TestWindow->setSizePolicy(sizePolicy);
        TestWindow->setMinimumSize(QSize(0, 0));
        TestWindow->setMaximumSize(QSize(16777215, 16777215));
        verticalLayout = new QVBoxLayout(TestWindow);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setSpacing(0);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        img1 = new QLabel(TestWindow);
        img1->setObjectName(QStringLiteral("img1"));

        horizontalLayout->addWidget(img1);

        img2 = new QLabel(TestWindow);
        img2->setObjectName(QStringLiteral("img2"));

        horizontalLayout->addWidget(img2);


        verticalLayout->addLayout(horizontalLayout);

        verticalLayout->setStretch(0, 8);

        retranslateUi(TestWindow);

        QMetaObject::connectSlotsByName(TestWindow);
    } // setupUi

    void retranslateUi(QDialog *TestWindow)
    {
        TestWindow->setWindowTitle(QApplication::translate("TestWindow", "Dialog", 0));
        img1->setText(QString());
        img2->setText(QString());
    } // retranslateUi

};

namespace Ui {
    class TestWindow: public Ui_TestWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_TESTWINDOW_H
