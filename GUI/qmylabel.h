// from qmylabel.py
/*
from PyQt4 import QtCore, QtGui, Qt
import default_parameters
##import ABM3
*/

//class QMyLabel(QtGui.QLabel):
//    '''
//    What I did was to set the labels as QMyLabel in QtDesigner, and reclass
//    its mousePressEvent.
//    '''
//    __pyqtSignals__ = ("labelClicked(Text)")

#ifndef QMYLABEL_H
#define QMYLABEL_H

#include <qlabel.h>
#include "log.h"
LOG_USE();

class QMyLabel: public QLabel 
{
    Q_OBJECT

public:
	QMyLabel(QWidget *parent = 0);
	
signals:
	void labelClicked(QString text);

private:

	void mousePressEvent (QMouseEvent *event);
	/*
	{
		QString text = "clicked";
//        emit(SIGNAL(labelClicked(QString)), text);
        emit(labelClicked(text));
		LOG_MSG("mousePressEvent emitted signal");
	};	
	*/
};

//QMyLabel::QMyLabel(QWidget *parent) : QLabel(parent)
//{}

#endif

    /*
    def mousePressEvent (self, event):
        """
##        defaultParameterList = [
##        ["TC_AVIDITY_MEDIAN", 1.0, 0.0, 10.0,
##        "TCR avidity median parameter",
##        "TCR avidity has a log normaldistribution, described by the median and shape parameters."],
##
         ...
##        ]
        

        paramName = self.objectName()[6:]
        for param in default_parameters.defaultParameterList: # We are only using the detailed info, not changable by the user, so using the defaultParamList will be fine.
            if param[0] == paramName:
                text = param[5]
                print "QMyLabel: ",text

##        ABM3.StartMain.showDescription()
        self.emit(QtCore.SIGNAL('labelClicked(QString)'), text)
		*/

//	void mousePressEvent (QEvent *event) {
//		QString text = "clicked";
//        emit(SIGNAL(labelClicked(QString)), text);
//	};
