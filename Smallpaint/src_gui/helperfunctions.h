#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H

#include <QString>
#include <QImage>
#include <main.h>

class HelperFunctions {

public:
	void saveImage();
	QString getRenderer(QString renderer);
	void changedComboBox(QString renderer);
	int getTextAsInt(QString text);
    float getTextAsFloat(QString text);
    void getInput(QString selectedRenderer, smallpaint::RenderInfo &info);
    void initializeRenderInformation(QString selectedRenderer, smallpaint::RenderInfo info);
	QString processTime(int time);
	void updateInfo(int currentSpp, int goalSpp);

private:

};

#endif // HELPERFUNCTIONS_H
