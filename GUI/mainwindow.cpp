/****************************************************************************
 ABM_GUI
****************************************************************************/

#include <QtGui>

#include "mainwindow.h"
#include "log.h"
#include "params.h"
#include "misc.h"
#include "plot.h"
#include "myvtk.h"
#include "transfer.h"
#ifdef _WIN32
#include "windows.h"
#define sleep(n) Sleep(1000 * n)
#endif

#ifdef __LINUX
#include <QTcpServer>
#else
#include <QTcpServer.h>
#endif

LOG_USE();

Params *parm;	// I don't believe this is the best way, but it works
int showingVTK;
int VTKbuffer[100];
int mono_list[5*MAX_MONO];
int nmono_list;
float cap_list[7*MAX_CAP];
int ncap_list;
float pit_list[4*MAX_PIT];
int npit_list;
float clast_list[7*MAX_CLAST];
int nclast_list;
int blast_list[5*MAX_BLAST];
int nblast_list;
QMutex mutex1, mutex2;

int summaryData[100];
int NX, NY, NZ, NBY;
int nt_vtk;
bool leftb;

QMyLabel::QMyLabel(QWidget *parent) : QLabel(parent)
{}

//--------------------------------------------------------------------------------------------------------
// Redefines mousePressEvent for QMyLabel, which extends QLabel.  This is used to display info about
// a model parameter.
//--------------------------------------------------------------------------------------------------------
void QMyLabel::mousePressEvent (QMouseEvent *event) {
	event->accept();
	QString sname = objectName().mid(6);
	QString text = "mousePressEvent";
	// Find which label_ sent the signal, and read its text
	for (int k=0; k<parm->nParams; k++) {
		PARAM_SET param = parm->get_param(k);
		if (sname.compare(param.tag) == 0)
			text = param.text;
	}
    emit labelClicked(text);
};	

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
MainWindow::MainWindow(QWidget *parent)
   : QMainWindow(parent)
{
	LOG_MSG("Started MainWindow");

	setupUi(this);
    showMaximized();
	// Some initializations
    nDistPts = 200;
	nTicks = 1000;
	tickVTK = 100;	// timer tick for VTK in milliseconds 
	paramSaved = false;
	paused = false;
	posdata = false;
    DCmotion = false;
    done = false;
    first = true;
	started = false;
    firstVTK = true;
	showingVTK = 0;
	nGraphCases = 0;
	for (int i=0; i<Plot::ncmax; i++) {
		graphResultSet[i] = 0;
	}
    cellfile = "cell_pos.dat";
	vtkfile = "basecase.pos";
	savepos_start = 0;
	ntimes = 0;
	hour = 0;

	DISABLE_TABS = true;

	QFileInfo cellfile_info(cellfile);
	LOG_QMSG(cellfile_info.path());
	LOG_QMSG(cellfile_info.fileName());
	if (QFile::exists(cellfile))
		QFile::remove(cellfile);
	param_to_sliderIndex = NULL;
	defaultInputFile = "basecase.inp";
	inputFile = defaultInputFile;

	parm = new Params();
	nParams = parm->nParams;
	createLists();
	createActions(); 
	loadParams();
	writeout();
    timer = new QTimer(this);
	vtk = new MyVTK(mdiArea_VTK);
	vtk->init();
	tabs->setCurrentIndex(0);
	goToInputs();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::createActions()
{
	action_stop->setEnabled(false);
    action_pause->setEnabled(false);
    action_inputs->setEnabled(false);
    action_outputs->setEnabled(false);
	action_save_snapshot->setEnabled(false);
	text_more->setEnabled(false);

    connect(action_open_input, SIGNAL(triggered()), this, SLOT(readInputFile()));
    connect(action_load_results, SIGNAL(triggered()), this, SLOT(loadResultFile()));
    connect(action_saveAs, SIGNAL(triggered()), this, SLOT(saveAs()));
    connect(action_save, SIGNAL(triggered()), this, SLOT(save()));
    connect(action_inputs, SIGNAL(triggered()), SLOT(goToInputs()));
    connect(action_outputs, SIGNAL(triggered()), SLOT(goToOutputs()));
    connect(action_VTK, SIGNAL(triggered()), SLOT(goToVTK()));
    connect(action_run, SIGNAL(triggered()), SLOT(runServer()));
    connect(action_pause, SIGNAL(triggered()), SLOT(pauseServer()));
    connect(action_stop, SIGNAL(triggered()), SLOT(stopServer()));
    connect(action_play_VTK, SIGNAL(triggered()), SLOT(playVTK()));
    connect(action_set_speed, SIGNAL(triggered()), SLOT(setVTKSpeed()));
	for (int i=0; i<nLabels; i++) {
		QLabel *label = label_list[i];
		QString label_str = label->objectName();
		if (label_str.startsWith("label_")) {
			connect((QObject *)label, SIGNAL(labelClicked(QString)), this, SLOT(showMore(QString)));
		}
	}
	// Graph menu
    connect(action_save_snapshot, SIGNAL(triggered()), this, SLOT(saveSnapshot()));
//	connect((QObject *)lcdNumber_hour, SIGNAL(hourUpdate(double)), this, SLOT(Display(double)));
	cbox_fastdisplay->setChecked(true);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::createLists()
{        
	lineEdit_list = findChildren<QLineEdit *>();
	spin_list = findChildren<QSpinBox *>();
	combo_list = findChildren<QComboBox *>();
	checkbox_list = findChildren<QCheckBox *>();
	radiobutton_list = findChildren<QRadioButton *>();
	slider_list = findChildren<QSlider *>();
	label_list = findChildren<QLabel *>();

	for (int i=0; i<lineEdit_list.length(); i++) {
		widget_list.append(lineEdit_list[i]);
	}
	for (int i=0; i<spin_list.length(); i++) {
		widget_list.append(spin_list[i]);
	}
	for (int i=0; i<combo_list.length(); i++) {
		widget_list.append(combo_list[i]);
	}
	for (int i=0; i<checkbox_list.length(); i++) {
		widget_list.append(checkbox_list[i]);
	}
	for (int i=0; i<radiobutton_list.length(); i++) {
		widget_list.append(radiobutton_list[i]);
	}

	nWidgets = widget_list.length();
	nSliders = slider_list.length();
	nLabels = label_list.length();
	sprintf(msg,"nLabels: %d",nLabels);
	LOG_MSG(msg);

	for (int i=0; i<nWidgets; i++) {
		QWidget *w = widget_list[i];
        QString wname = w->objectName();
		if (wname.startsWith("line_")) {
			connect(w, SIGNAL(textChanged(QString)), this, SLOT(changeParam()));
		}
		if (wname.startsWith("spin_")) {
			connect(w, SIGNAL(valueChanged(int)), this, SLOT(changeParam()));
		}
		if (wname.startsWith("comb_")) {
			connect(w, SIGNAL(activated(QString)), this, SLOT(changeParam()));
		}
		if (wname.startsWith("cbox_")) {
			connect(w, SIGNAL(toggled(bool)), this, SLOT(changeParam()));
		}
		if (wname.startsWith("rbut_")) {
			connect(w, SIGNAL(toggled(bool)), this, SLOT(changeParam()));
		}
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: drawDistPlots()
{
	double *x, *prob;
	x = new double[nDistPts];
	prob = new double[nDistPts];
	QwtPlot *qp;
	string name_str;
	QString median_qstr, shape_qstr;
	double median, shape;
	for (int j=0; j<3; j++) {
		qp = distplot_list[j];
        QString name = qp->objectName();
		LOG_QMSG(name);
		LOG_MSG("drawDistPlots (a)")
		median = median_qstr.toDouble();
		shape = shape_qstr.toDouble();
        create_lognorm_dist(median,shape,nDistPts,x,prob);
		LOG_MSG("drawDistPlots (b)")

        int n = dist_limit(prob,nDistPts);
        double xmax = x[n];
		sprintf(msg,"%f %f %d",median,shape,n);
		for (int i=0;i<40;i++) {
			sprintf(msg,"%d %f %f",i,x[i],prob[i]);
		}
        qp->setAxisScale(QwtPlot::xBottom, 0.0, xmax, 0.0);
        QwtPlotCurve *curve = new QwtPlotCurve("title");
		LOG_MSG("drawDistPlots (c)")
		curve->attach(qp);
        curve->setData(x, prob, n);
        curve_list[j] = curve;
		qp->replot();
		LOG_MSG("drawDistPlots (d)")
	}
	delete [] x;
	x = NULL;
	delete [] prob;
	prob = NULL;
}

//-------------------------------------------------------------
// Loops through the workingParameterList and fills in the GUI.
//-------------------------------------------------------------
void MainWindow::loadParams()
{
	sprintf(msg,"nParams: %d nSliders: %d",nParams,nSliders);
	for (int k=0; k<nSliders; k++) {		// there must be a neater way to do this
		SliderPlus *splus = 0;
		sliderplus_list.append(splus);
		QWidget *w = 0;
		sliderParam.append(w);		
	}
	if (param_to_sliderIndex == NULL) {
		param_to_sliderIndex = new int[nParams];
		for (int i=0; i<nParams; i++)
			param_to_sliderIndex[i] = -1;
	}
	for (int i=0; i<nWidgets; i++) {
		QWidget *w = widget_list[i];							// w = widget_list[i] is the ith widget in the UI
        QString qsname = w->objectName();
		if (qsname.startsWith("line_") || qsname.startsWith("spin_") 
			|| qsname.startsWith("comb_") || qsname.startsWith("cbox_")
			|| qsname.startsWith("rbut_")) {
			QString wtag = qsname.mid(5);
			int rbutton_case = 0;
			if (qsname.startsWith("rbut_")) {
				wtag = parse_rbutton(wtag,&rbutton_case);
			}
            // Find corresponding data in workingParameterList
            bool found = false;
			for (int k=0; k<nParams; k++) {
				PARAM_SET p = parm->get_param(k);
				QString ptag = p.tag;		// ptag is the tag of the kth parameter in the list
				if (wtag.compare(ptag) == 0) {
					double vmax = p.maxvalue;
					double vmin = p.minvalue;
                    // Update the widget (line_, spin_ or comb_) with data from the parameter list
                    // ---LineEdits
					if (qsname.startsWith("line_")) {
                        double val = p.value;
						QString val_str = QString::number(val);
						QLineEdit *w_l = (QLineEdit *)w;
                        w_l->setText(val_str);
						if (USE_RANGES) {
							// Set max and min values. If min=max=0, there're no restrictions.
							if (!(vmin == 0 && vmax == 0)) {
								QValidator *aValidator = new MyDoubleValidator(vmin, vmax, 8, w_l);
								w_l->setValidator(aValidator);
							}
						}
					} else if (qsname.startsWith("spin_")) {
                        double val = p.value;
						QSpinBox *w_s = (QSpinBox *)w;
                        w_s->setValue(val);
						if (!(vmin == 0 && vmax == 0)) {
                            w_s->setMinimum(vmin);
                            w_s->setMaximum(vmax);
						}
						if (qsname.contains("NCPU")) {
							ncpu = p.value;
						}
					} else if (qsname.startsWith("comb_")) {
                        int val = p.value - 1;	//0-based indexing
						QComboBox *w_c = (QComboBox *)w;
                        w_c->setCurrentIndex(val);
                        /*
					} else if (qsname.startsWith("cbox_")) {
						QCheckBox *w_cb = (QCheckBox *)w;
						bool in_vitro = qsname.contains("IN_VITRO");
						if (p.value == 1) {
							w_cb->setChecked(true);
							if (in_vitro) 
								enableInVitro();
						} else {
							w_cb->setChecked(false);
							if (in_vitro) 
								disableInVitro();
						}
						bool dc_injection = qsname.contains("DC_INJECTION");
						if (p.value == 1) {
							w_cb->setChecked(true);
							if (dc_injection)
								enableDCInjection();
						} else {
							w_cb->setChecked(false);
							if (dc_injection)
								disableDCInjection();
						}
						bool use_traffic = qsname.contains("USE_TRAFFIC");
						if (p.value == 1) {
							w_cb->setChecked(true);
							if (use_traffic)
								enableUseTraffic();
						} else {
							w_cb->setChecked(false);
							if (use_traffic)
								disableUseTraffic();
						}
                        */
					} else if (qsname.startsWith("rbut_")) {
						QRadioButton *w_rb = (QRadioButton *)w;
						if (p.value == rbutton_case) {
							w_rb->setChecked(true);
						} else {
							w_rb->setChecked(false);
						}
					}
					
                    // Update Label text
                    // Get the corresponding label from the label list
                    QString labelString = "label_" + wtag;
                                        QLabel *label = NULL;
					bool foundLabel = false;
					for (int j=0; j<nLabels; j++) {
						label = label_list[j];
						if (labelString.compare(label->objectName()) == 0) {
							foundLabel = true;
							break;
						}
					}										// label is the pointer to the UI label for wtag and ptag
                    QString labelText = p.label;
                    
                    // Hardcode the distribution label names for now
                    if (wtag.compare("TC_AVIDITY_MEDIAN") == 0)
                        labelText = "Median";
					else if (wtag.compare("TC_AVIDITY_SHAPE") == 0)
                        labelText = "Shape";
                    else if (wtag.compare("DIVIDE1_MEDIAN") == 0)
                        labelText = "Median";
					else if (wtag.compare("DIVIDE1_SHAPE") == 0)
                        labelText = "Shape";
                    else if (wtag.compare("DIVIDE2_MEDIAN") == 0)
                        labelText = "Median";
					else if (wtag.compare("DIVIDE2_SHAPE") == 0)
                        labelText = "Shape";
                    else if (wtag.compare("DC_ANTIGEN_MEDIAN") == 0)
                        labelText = "Median";
					else if (wtag.compare("DC_ANTIGEN_SHAPE") == 0)
                        labelText = "Shape";
                    else if (wtag.compare("DC_LIFETIME_MEDIAN") == 0)
                        labelText = "Median";
					else if (wtag.compare("DC_LIFETIME_SHAPE") == 0)
                        labelText = "Shape";

					bool is_slider = false;
					int j;
					QSlider *s;
					QString sliderString;
					for (j=0; j<nSliders; j++) {
						sliderString = "slider_" + wtag;
						s = slider_list[j];
						if (sliderString.compare(s->objectName()) == 0) {
							is_slider = true;					// the jth slider in the list corresponds to wtag and ptag
							break;
						}
					}
					if (is_slider) {
                        // If there is a slider corresponding to wtag, then just use the label.
						if (foundLabel)
	                        label->setText(labelText);
					} else {
						if (!(vmin == 0 && vmax == 0)) {
                            // If there is no slider, then add min and max values to the label text.
							QString min_str = QString::number(vmin);
							QString max_str = QString::number(vmax);
							if (foundLabel)
		                        label->setText(labelText + "  [ " + min_str + "-" + max_str + " ]");
						} else {
							if (foundLabel)
		                        label->setText(labelText);
						}
					}
						
                    // If there is a corresponding slider for this parameter, then apply settings.
					if (is_slider) {						
                        SliderPlus *splus = new SliderPlus(wtag,vmin,vmax,nTicks,k,i);
						sliderplus_list[j] = splus;
                        int ival = splus->val_to_int(p.value);
                        s->setMinimum(0);
                        s->setMaximum(splus->nTicks());
                        s->setSliderPosition(ival);
						sliderParam[j] = w;
                        connect(s, SIGNAL(valueChanged(int)), this, SLOT(updateSliderBox())); //sliderReleased
                        
                        param_to_sliderIndex[k] = j;
					}                  
                    found = true;
                    break;
				}
			}
		}
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
QString MainWindow::parse_rbutton(QString wtag, int *rbutton_case)
{
	// parse wtag into part before '_' and part after '_'
	int j = wtag.indexOf('_');
	QString suffix = wtag.mid(j+1);
	// the prefix becomes wtag, the suffix becomes rbutton_case, an integer 0,1,2,...
	wtag = wtag.mid(0,j);
	bool ok;
	*rbutton_case = suffix.toInt(&ok);
	return wtag;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::reloadParams()
{
	for (int i=0; i<nWidgets; i++) {
		QWidget *w = widget_list[i];							// w = widget_list[i] is the ith widget in the UI
        QString qsname = w->objectName();
		if (qsname.startsWith("line_") || qsname.startsWith("spin_") 
			|| qsname.startsWith("comb_") || qsname.startsWith("cbox_")
			|| qsname.startsWith("rbut_")) {
//			LOG_QMSG(qsname);
			QString wtag = qsname.mid(5);
			int rbutton_case = 0;
			if (qsname.startsWith("rbut_")) {
				wtag = parse_rbutton(wtag,&rbutton_case);
			}
            // Find corresponding data in workingParameterList
            bool found = false;
			for (int k=0; k<nParams; k++) {
				PARAM_SET p = parm->get_param(k);
				QString ptag = p.tag;		// ptag is the tag of the kth parameter in the list
				if (wtag.compare(ptag) == 0) {
					found = true;
                    // Update the widget (line_, spin_ or comb_) with data from the parameter list
					if (qsname.startsWith("line_")) {
                        double val = p.value;
						QString val_str = QString::number(val);
						QLineEdit *w_l = (QLineEdit *)w;
                        w_l->setText(val_str);
					} else if (qsname.startsWith("spin_")) {
                        double val = p.value;
						QSpinBox *w_s = (QSpinBox *)w;
                        w_s->setValue(val);
						if (qsname.contains("NCPU")) {
							ncpu = p.value;
						}
					} else if (qsname.startsWith("comb_")) {
                        int val = p.value - 1;	//0-based indexing
						QComboBox *w_c = (QComboBox *)w;
                        w_c->setCurrentIndex(val);
                        /*
					} else if (qsname.startsWith("cbox_")) {
						QCheckBox *w_cb = (QCheckBox *)w;
						bool in_vitro = qsname.contains("IN_VITRO");
						if (p.value == 1) {
							w_cb->setChecked(true);
							if (in_vitro) 
								enableInVitro();
						} else {
							w_cb->setChecked(false);
							if (in_vitro) 
								disableInVitro();
						}
						bool dc_injection = qsname.contains("DC_INJECTION");
						if (p.value == 1) {
							w_cb->setChecked(true);
							if (dc_injection)
								enableDCInjection();
						} else {
							w_cb->setChecked(false);
							if (dc_injection)
								disableDCInjection();
						}
						bool use_traffic = qsname.contains("USE_TRAFFIC");
						if (p.value == 1) {
							w_cb->setChecked(true);
							if (use_traffic)
								enableUseTraffic();
						} else {
							w_cb->setChecked(false);
							if (use_traffic)
								disableUseTraffic();
						}
                        */
					} else if (qsname.startsWith("rbut_")) {
						QRadioButton *w_rb = (QRadioButton *)w;
						if (p.value == rbutton_case) {
							w_rb->setChecked(true);
						} else {
							w_rb->setChecked(false);
						}
					}
				}
			}
			if (!found) {
				LOG_MSG("Widget tag not found:");
				LOG_QMSG(qsname);
				LOG_QMSG(wtag);
			}
		}
	}				
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::showMore(QString moreText)
{
	LOG_MSG("label clicked!");
	LOG_QMSG(moreText);
	
//	if ((int)sender() != currentDescription) {
	if ((uintptr_t)sender() != currentDescription) {
		text_more->setEnabled(true); // self.ui.text_description.setEnabled(1) #show()
        text_more->setText(moreText); // text_description
		currentDescription = (uintptr_t)sender();
    } else {
        text_more->clear(); // text_description
        text_more->setEnabled(false); // hide()#text_description
        currentDescription = 0;
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::writeout()
{
	QString line;
    QFile file(inputFile);
	LOG_MSG("writeout");
	LOG_QMSG(inputFile);
    if (!file.open(QFile::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot write file %1:\n%2.")
                             .arg(inputFile)
                             .arg(file.errorString()));
		LOG_MSG("File open failed");
        return;
    }
    QTextStream out(&file);
	for (int k=0; k<parm->nParams; k++) {
		PARAM_SET p = parm->get_param(k);
		double val = p.value;
		if (val == int(val)) 	// whole number, write as integer
			line = QString::number(int(val));
		else
			line = QString::number(val);
		int nch = line.length();
		for (int i=0; i<max(12-nch,1); i++)
			line += " ";
		line += p.tag;
		line += "\n";
		out << line;
	}

    paramSaved = true;
	LOG_MSG("Input data saved");
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::readInputFile()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open ..."), ".", tr("Input Files (*.inp)"));
	if (fileName.compare("") == 0)
		return;
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot read file %1:\n%2.")
                             .arg(fileName)
                             .arg(file.errorString()));
        return;
    }

    QTextStream in(&file);
	QString line;
	for (int k=0; k<parm->nParams; k++) {
		line = in.readLine();
		QStringList data = line.split(" ",QString::SkipEmptyParts);
		parm->set_value(k,data[0].toDouble());
	}

    reloadParams();
    paramSaved = true;
	inputFile = fileName;
}

//--------------------------------------------------------------------------------------------------------
// NOT USED
//--------------------------------------------------------------------------------------------------------
void MainWindow::loadResultFile()
{
	RESULT_SET *R;

	R = new RESULT_SET;
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open ..."), ".", tr("Result Files (*.res)"));
	if (fileName.compare("") == 0)
		return;
	R->casename = QFileInfo(fileName).baseName();
	for(int i=0; i<result_list.size(); i++) {
		if (R->casename.compare(result_list[i]->casename) == 0) {
			QMessageBox::warning(this, tr("Open results"),
                             tr("This result file is already loaded"));
			return;
		}
	}
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot read file %1:\n%2.")
                             .arg(fileName)
                             .arg(file.errorString()));
        return;
    }

    QTextStream in(&file);
	R->nsteps = 0;
	bool indata = false;
	QString line;
	do {
		line = in.readLine();
		if (line.length() > 0) {
			QStringList datalist = line.split(" ",QString::SkipEmptyParts);
			if (indata) {
				R->nsteps++;
			}
			if (datalist[0].contains("nlist")) {
				indata = true;
			}
		}
	} while (!line.isNull());
	in.seek(0);
	R->tnow = new double[R->nsteps];
	R->nDC = new double[R->nsteps];
	R->nborn = new double[R->nsteps];
	R->ntot = new double[R->nsteps];
	R->ncogseed = new double[R->nsteps];
	R->ncog = new double[R->nsteps];
	R->ndead = new double[R->nsteps];
	R->teffgen = new double[R->nsteps];

	step = -1;
	indata = false;
	do {
		line = in.readLine();
		if (line.length() > 0) {
			QStringList dataList = line.split(" ",QString::SkipEmptyParts);
			if (indata) {
				double data[10];
				for (int k=0; k<10; k++) 
					data[k] = dataList[k].toDouble();
				step++;
				if (step >= R->nsteps) {
					LOG_MSG("ERROR: loadResultFile: step >= nsteps_p");
					return;
				}
				R->tnow[step] = step;		//data[1];
				R->nDC[step] = data[2];
				R->nborn[step] = data[3];
				R->ntot[step] = data[4];
				R->ncogseed[step] = data[5];
				R->ncog[step] = data[6];
				R->ndead[step] = data[8];
				R->teffgen[step] = data[9];
			}
			if (dataList[0].contains("nlist")) {
				indata = true;
			}
		}
	} while (!line.isNull());

	// Compute the maxima
	R->max_nborn = getMaximum(R,R->nborn);
	R->max_ncog = getMaximum(R,R->ncog);
	R->max_ncogseed = getMaximum(R,R->ncogseed);
	R->max_nDC = getMaximum(R,R->nDC);
	R->max_teffgen = getMaximum(R,R->teffgen);
	R->max_ntot = getMaximum(R,R->ntot);

	// Now add the result set to the list
	result_list.append(R);
	if (nGraphCases == 0) {
		if (show_outputdata)
			box_outputData = 0;
		initializeGraphs(R);
		drawGraphs();
		goToOutputs();
		sprintf(msg,"nGraphCases: %d %p",nGraphCases,graphResultSet[0]);
		LOG_MSG(msg);
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
bool MainWindow::save()
{
	writeout();
	return true;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
bool MainWindow::saveAs()
{
    // show the file dialog
	QString fileName = QFileDialog::getSaveFileName(this, tr("Select Input File"), ".", tr("Input Files (*.inp)"));    
	if (fileName.compare("") != 0) {
		LOG_MSG("Selected file:");
		LOG_QMSG(fileName);
		inputFile = fileName;
        writeout();
	}
    // Otherwise if user chooses cancel ...
	return true;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
double MainWindow::getMaximum(RESULT_SET *R, double *x)
{
	double maxx = 0;
	for (int i=0; i<R->nsteps; i++)
		maxx = max(maxx,x[i]);
	return maxx;
}

//-------------------------------------------------------------
// Switches to the input screen
// For when the outputs are being displayed
//-------------------------------------------------------------
void MainWindow::goToInputs()
{
    stackedWidget->setCurrentIndex(0);
	showingVTK = 0;
	action_inputs->setEnabled(false);
    action_outputs->setEnabled(true);
    action_VTK->setEnabled(true);
}

//-------------------------------------------------------------
// Switches to the output screen
//-------------------------------------------------------------
void MainWindow::goToOutputs()
{
    stackedWidget->setCurrentIndex(1);    
	showingVTK = 0;
	action_outputs->setEnabled(false);
    action_inputs->setEnabled(true);
    action_VTK->setEnabled(true);
}

//-------------------------------------------------------------
// Switches to the VTK screen
//-------------------------------------------------------------
void MainWindow::goToVTK()
{
	LOG_QMSG("goToVTK");
	if (started) {
		stackedWidget->setCurrentIndex(2);
		action_outputs->setEnabled(true);
		action_inputs->setEnabled(true);
		action_VTK->setEnabled(false);
		showingVTK = 1;
		LOG_QMSG("showingVTK");
	}
}

//-------------------------------------------------------------
// Load and play stored cell position data
//-------------------------------------------------------------
void MainWindow::playVTK()
{
	LOG_MSG("playVTK");
	// Select a file to play
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open ..."), ".", tr("Cell Path Files (*.pos)"));
	if (fileName.compare("") == 0)
		return;
    QMessageBox::StandardButton reply;
    reply = QMessageBox::question(this, tr("Animation player"),
									tr("Save animation frames to image files (.jpg)?"),
                                    QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel);
	bool save_image;
	if (reply == QMessageBox::Yes)
		save_image = true;
	else
		save_image = false;
	started = true;
	showingVTK = 0;
	goToVTK();
	if (!vtk->startPlayer(QFileInfo(fileName).absoluteFilePath(), timer, save_image)) {
		LOG_MSG("startPlayer failed");
		errorPopup("Open failure on this file");
		return;
//		exit(1);
	}
    connect(timer, SIGNAL(timeout()), this, SLOT(timer_update()));
    timer->start(tickVTK);
	action_stop->setEnabled(true);
    action_pause->setEnabled(true);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::timer_update()
{
	ntimes ++;
	if (!vtk->nextFrame()) {
		LOG_MSG("Player completed");
		action_stop->setEnabled(false);
		action_pause->setEnabled(false);
//		timer->stop();
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::setVTKSpeed()
{
    bool ok;
    int i = QInputDialog::getInt(this, tr("Set speed"),
		tr("Player timer tick (ms): "), tickVTK, 10, 10000, 1, &ok);
	if (ok) {
		tickVTK = i;
	}
	if (started) {
		timer->stop();
		timer->start(tickVTK);
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::setSavePosStart()
{
    bool ok;
    int i = QInputDialog::getInt(this, tr("Set savepos start"),
		tr("Start recording cell positions at (hours): "), savepos_start, 0, 1000, 1, &ok);
	if (ok) {
		savepos_start = i;
		sprintf(msg,"savepos_start: %d",savepos_start);
		LOG_MSG(msg);
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::saveSnapshot()
{
	LOG_MSG("saveSnapshot");
	QString fileName = QFileDialog::getSaveFileName(this, tr("Select image file"), ".", 
		tr("Image files (*.png *.jpg *.tif *.bmp)"));    
	if (fileName.compare("") == 0) {
		goToVTK();
		return;
	}
	QFileInfo fi(fileName);
	QString imgType = fi.suffix();
	goToVTK();
	vtk->saveSnapshot(fileName,imgType);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::runServer()
{
	if (paused) {
		if (vtk->playing) {
			vtk->playon();
		} else {
			exthread->unpause();
		}
        action_run->setEnabled(false);
        action_pause->setEnabled(true);
        action_stop->setEnabled(true);
		action_save_snapshot->setEnabled(false);
        paused = false;
        return;
	}
	if (!paramSaved) {
		int response = QMessageBox::critical(this, tr("ABM Model GUI"), \
					tr("The document has been modified.\nPlease save changes before continuing."), \
					QMessageBox::Save | QMessageBox::Cancel); // | Qt.QMessageBox.Discard
		if (response == QMessageBox::Save) {
            save();
		} else if (response == QMessageBox::Cancel) {
            return;
		}
	}

	if (!first) {
		int response = QMessageBox::question(this, tr("ABM Model GUI"),
						tr("Would you like to clear the graphs from the previous run?"),
						QMessageBox::Yes | QMessageBox::No);
		if (response == QMessageBox::Yes)
			mdiArea_plots->closeAllSubWindows();
		else if (response == QMessageBox::Cancel)
            return;
	}

    // Display the outputs screen
	if (stackedWidget->currentIndex() == 2) {
		showingVTK = 1;
	} else {
		stackedWidget->setCurrentIndex(1);
		showingVTK = 0;
	}

    // Disable parts of the GUI        
    action_run->setEnabled(false);
	action_pause->setEnabled(false);
	action_stop->setEnabled(false);
    action_inputs->setEnabled(true);
    action_VTK->setEnabled(true);
	action_save_snapshot->setEnabled(false);
	tab_monocyte->setEnabled(false);
	tab_osteoclast->setEnabled(false);
	tab_signal->setEnabled(false);
    tab_run->setEnabled(false);

	if (show_outputdata)
	    box_outputData = new QTextBrowser();
	else
		box_outputData = 0;

	if (use_CPORT1) {

		if (QFile::exists(cellfile))
			QFile::remove(cellfile);

		// Port 5001
		sthread1 = new SocketHandler(CPORT1);
		connect(sthread1, SIGNAL(sh_output(QString)), this, SLOT(outputData(QString)));
		sthread1->start();
	}

	// Port 5000
	sthread0 = new SocketHandler(CPORT0);
	connect(sthread0, SIGNAL(sh_output(QString)), box_outputLog, SLOT(append(QString))); //self.outputLog)
	connect(sthread0, SIGNAL(sh_connected()), this, SLOT(preConnection()));
	connect(sthread0, SIGNAL(sh_disconnected()), this, SLOT(postConnection()));
	sthread0->start();
	vtk->cleanup();
	Sleep(100);

	hours = 0;
	nt_vtk = 0;
	for (int k=0; k<parm->nParams; k++) {
		PARAM_SET p = parm->get_param(k);
		if (p.tag.compare("NDAYS") == 0) {
			hours = p.value*24;
		}
		if (p.tag.compare("NT_ANIMATION") == 0) {
			nt_vtk = p.value;
		}
	}

    started = true;
    if (rbut_RUNCASE_0->isChecked())
        runcase = MONO_CASE;
    else
        runcase = OC_CASE;
    exthread = new ExecThread(runcase,inputFile);
	connect(exthread, SIGNAL(display()), this, SLOT(displayScene()));
	connect(exthread, SIGNAL(summary()), this, SLOT(showSummary()));
	connect(exthread, SIGNAL(initialized()), this, SLOT(setInitialized()));
	exthread->ncpu = ncpu;
	exthread->nsteps = int(hours*60/DELTA_T);
	exthread->paused = false;
	exthread->stopped = false;
	exthread->start();

	goToVTK();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::setInitialized()
{
	action_pause->setEnabled(true);
	action_stop->setEnabled(true);
}

void MainWindow::zzz()
{
	LOG_MSG("zzz");
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::preConnection()
{
	LOG_MSG("preConnection");
	// We assume that the model output is at hourly intervals
	newR = new RESULT_SET;
	QString casename = QFileInfo(inputFile).baseName();
	vtkfile = casename + ".pos";
	newR->casename = casename;
	int nsteps = int(hours+1.5);
	sprintf(msg,"preConnection: nsteps: %d",nsteps);
	LOG_MSG(msg);
	newR->nsteps = nsteps;
	newR->tnow = new double[nsteps];
	newR->nDC = new double[nsteps];
	newR->nborn = new double[nsteps];
	newR->ntot = new double[nsteps];
	newR->ncogseed = new double[nsteps];
	newR->ncog = new double[nsteps];
	newR->ndead = new double[nsteps];
	newR->teffgen = new double[nsteps];
	LOG_MSG("preconnection: Allocated result set arrays");

	newR->tnow[0] = 0;	// These are not the right initial values
	newR->nDC[0] = 0;
	newR->nborn[0] = 0;
	newR->ntot[0] = 0;
	newR->ncogseed[0] = 0;
	newR->ncog[0] = 0;
	newR->ndead[0] = 0;
	newR->teffgen[0] = 0;

    step = -1;
    // Initialize graphs
	initializeGraphs(newR);
    posdata = false;
	if (cbox_savepos->isChecked()) {
		if (QFile::exists(vtkfile))
			QFile::remove(vtkfile);
	}
	LOG_MSG("preconnection: done");
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::errorPopup(QString errmsg)
{
	LOG_QMSG(errmsg);
	QMessageBox msgBox;
	msgBox.setText(errmsg);
	msgBox.exec();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::initializeGraphs(RESULT_SET *R)
{
	mdiArea_plots->closeAllSubWindows();
	mdiArea_plots->show();
	clearAllGraphs();

	graph_ntot = new Plot("ntot",R->casename);
	graph_ntot->setTitle("Total Monocyte Population");
	graph_ntot->setAxisTitle(QwtPlot::yLeft, "No. of Cells");
	mdiArea_plots->addSubWindow(graph_ntot);
	graph_ntot->show();
	graph_ntot->setAxisScale(QwtPlot::xBottom, 0, R->nsteps, 0);

	graph_nborn = new Plot("nborn",R->casename);
	graph_nborn->setTitle("Total Monocytes Created");
	graph_nborn->setAxisTitle(QwtPlot::yLeft, "No. of Cells");
	mdiArea_plots->addSubWindow(graph_nborn);
	graph_nborn->show();
	graph_nborn->setAxisScale(QwtPlot::xBottom, 0, R->nsteps, 0);

	mdiArea_plots->tileSubWindows();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::drawGraphs()
{
	RESULT_SET *R;
	double nborn_max = 0, ntot_max = 0;
	for (int k=0; k<Plot::ncmax; k++) {
		R = graphResultSet[k];
		if (R != 0) {
			graph_ntot->redraw(R->tnow, R->ntot, R->nsteps, R->casename);
			ntot_max = max(ntot_max,R->max_ntot);
			nborn_max = max(nborn_max,R->max_nborn);
		}
	}
	graph_ntot->setYScale(ntot_max);
	graph_ntot->replot();
	graph_nborn->setYScale(nborn_max);
	graph_nborn->replot();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::displayScene()
{
	started = true;
	mutex2.lock();
	bool fast = cbox_fastdisplay->isChecked();
	vtk->get_cell_positions(fast);
	vtk->renderCells();
	mutex2.unlock();

}


//--------------------------------------------------------------------------------------------------------
// Currently summaryData[] holds istep,ntot,nborn.  Hourly intervals, i.e. every 240 timesteps
//--------------------------------------------------------------------------------------------------------
void MainWindow::showSummary()
{
//	LOG_MSG("showSummary");
	step++;
	if (step >= newR->nsteps) {
		LOG_MSG("ERROR: step >= nsteps");
		return;
	}
	mutex1.lock();
	hour = summaryData[0]*DELTA_T/60;
	progress = int(100.*hour/hours);
	progressBar->setValue(progress);

	QString hourstr = QString::number(int(hour));
	hour_display->setText(hourstr);

	QString casename = newR->casename;
	newR->tnow[step] = step;
	newR->ntot[step] = summaryData[1];
	newR->nborn[step] = summaryData[2];
	mutex1.unlock();
	graph_ntot->redraw(newR->tnow, newR->ntot, step+1, casename);
	graph_nborn->redraw(newR->tnow, newR->nborn, step+1, casename);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::outputData(QString qdata)
{
	if (qdata.startsWith("VTK")) {
		qdata.replace(0,3,"");
		bool savepos = cbox_savepos->isChecked();
		if (savepos) {
			if (step < savepos_start) {
				savepos = false;
			}
		}
		vtk->read_cell_positions(cellfile, vtkfile, savepos);
		started = true;
		if (showingVTK > 0 || firstVTK) {
			firstVTK = false;
			bool redo = false;
			if (showingVTK == 1) {
				redo = true;
				showingVTK = 2;
			} 
//			vtk->renderCells(redo,false);
		}
	    posdata = true;
		if (qdata.length() == 0)
			return;
	}
	if (quitMessage(qdata) || qdata.contains("Fortran") ) {
		return;
	}
	if (show_outputdata)
	    box_outputData->append(qdata);

    QStringList dataList = qdata.split(" ",QString::SkipEmptyParts);
	int data[4];
	for (int k=0; k<4; k++)
		data[k] = dataList[k].toInt();
	step++;
	if (step >= newR->nsteps) {
		LOG_MSG("ERROR: step >= nsteps");
		return;
	}
	hour = data[0]*DELTA_T/60;
	progress = int(100.*hour/hours);
	progressBar->setValue(progress);

// lcdNumber didn't work very well - too faint
//	lcdNumber_hour->display(int(hour));
//	lcdNumber_hour->show();
	QString hourstr = QString::number(int(hour));
	hour_display->setText(hourstr);

	QString casename = newR->casename;
    newR->tnow[step] = step;		//data[1];
	newR->ntot[step] = data[1];
	graph_ntot->redraw(newR->tnow, newR->ntot, step+1, casename);
	newR->nborn[step] = data[2];
	graph_nborn->redraw(newR->tnow, newR->nborn, step+1, casename);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::postConnection()
{
	LOG_MSG("postConnection");
	if (use_CPORT1) {
		sthread1->socket->close();
		sthread1->tcpServer->close();
		sthread1->quit();
		sthread1->wait(100);
		if (sthread1->isRunning()) {
			LOG_MSG("sthread1 did not terminate");
		}
	}

    action_run->setEnabled(true);
    action_pause->setEnabled(false);
    action_stop->setEnabled(false);
	action_save_snapshot->setEnabled(true);
	if (!DISABLE_TABS) {
	}
	tab_monocyte->setEnabled(true);
	tab_osteoclast->setEnabled(true);
	tab_signal->setEnabled(true);
	tab_run->setEnabled(true);

	// Check if a result set of this name is already in the list, if so remove it
	for (int i=0; i<result_list.size(); i++) {
		if (newR->casename.compare(result_list[i]->casename) == 0) {
			result_list.removeAt(i);
		}
	}
    posdata = false;
	LOG_MSG("completed postConnection");
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::close_sockets()
{
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::pauseServer()
{
	if (vtk->playing) {
		vtk->pause();
		LOG_MSG("Paused the player");
	} else {
		exthread->pause();
		LOG_MSG("Paused the ABM program.");
	}
	paused = true;
	action_run->setEnabled(true); 
	action_pause->setEnabled(false);
	action_stop->setEnabled(true);
	action_save_snapshot->setEnabled(true);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::stopServer()
{
	if (vtk->playing) {
		vtk->stop();
		LOG_MSG("Stopped the player");
	} else {
		LOG_MSG("stop ordered");
		if (paused) {
			LOG_MSG("was paused, runServer before stopping");
			runServer();
		}
		exthread->snapshot();
		LOG_MSG("exthread->stop");
		exthread->stop();
		LOG_MSG("exthread->stop return");
		sleep(1);		// delay for Fortran to wrap up (does this help?)
		if (use_CPORT1) {
			sthread1->quit();
			sthread1->terminate();
		}
		sthread0->stop();
		newR->nsteps = step+1;
	}
    action_run->setEnabled(true); 
    action_pause->setEnabled(false);
	action_stop->setEnabled(false);
	action_save_snapshot->setEnabled(true);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::clearAllGraphs()
{
	if (nGraphCases > 0) {
		graph_nborn->removeAllCurves();
		graph_ntot->removeAllCurves();
		graph_nDC->removeAllCurves();
		graph_teffgen->removeAllCurves();
		graph_ncog->removeAllCurves();
		if (!show_outputdata)
			graph_ncogseed->removeAllCurves();
		nGraphCases = 0;
	}
	for (int i=0; i<Plot::ncmax; i++) {
		graphCaseName[i] = "";
		graphResultSet[i] = 0;
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
QString MainWindow::selectResultSet()
{
    QStringList items;
	for (int k=0; k<result_list.size(); k++) {
		RESULT_SET *R = result_list[k];
		if (R == 0) continue;
		bool inlist = false;
		for (int i=0; i<Plot::ncmax; i++) {
			if (graphResultSet[i] == 0) continue;
			if (R->casename.compare(graphResultSet[i]->casename) == 0) {
				inlist = true;
				break;
			}
		}
		if (!inlist)
			items << R->casename;
	}
	if (items.size() == 0) {
		QMessageBox::warning(this, tr("Select result case"),
			tr("No result sets available - use 'File > Load results'"));
		return QString("");
	}

    bool ok;
    QString item = QInputDialog::getItem(this, tr("Select result case"),
		tr("Case:"), items, 0, false, &ok);
    if (ok && !item.isEmpty())
		return item;
	else
		return QString("");
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::addGraph()
{
	// Need to select a result set from result_list then add the corresponding curves
        RESULT_SET *R = NULL;

	if (nGraphCases == Plot::ncmax) {
		QString mess = QString("The maximum number of cases is %1").arg(Plot::ncmax);
		QMessageBox::warning(this, tr("Add graph"),	tr((mess.toStdString()).data())); 
		return;
	}
	QString casename = selectResultSet();
	if (casename.compare("") == 0)
		return;

	for (int k=0; k<result_list.size(); k++) {
		if (casename.compare(result_list[k]->casename) == 0) {
            R = result_list[k];		// OK after doing a run or a load, followed by another load
			break;
		}
	}

	graphResultSet[nGraphCases] = R;
	nGraphCases++;
	// First add the curves
	graph_nborn->addCurve(R->casename);
	graph_ntot->addCurve(R->casename);
	graph_nDC->addCurve(R->casename);
	graph_teffgen->addCurve(R->casename);
	graph_ncog->addCurve(R->casename);
	if (!show_outputdata)
		graph_ncogseed->addCurve(R->casename);
	// Adjust the x axis scale
	graph_nborn->setAxisAutoScale(QwtPlot::xBottom);
    graph_ntot->setAxisAutoScale(QwtPlot::xBottom);
    graph_ncog->setAxisAutoScale(QwtPlot::xBottom);
    graph_nDC->setAxisAutoScale(QwtPlot::xBottom);
    graph_teffgen->setAxisAutoScale(QwtPlot::xBottom);
	if (!show_outputdata)
	    graph_ncogseed->setAxisAutoScale(QwtPlot::xBottom);
	// Now redraw with the data
	drawGraphs();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
int MainWindow::selectGraphCase()
{
    QStringList items;
	for (int i=0; i<Plot::ncmax; i++) {
		if (graphResultSet[i] == 0) continue;
		items << graphResultSet[i]->casename;
	}
	if (items.size() == 0) {
		QMessageBox::warning(this, tr("Select graph case"),
			tr("No graph cases to remove"));
		return -1;
	}

    bool ok;
    QString item = QInputDialog::getItem(this, tr("Select graph case"),
		tr("Case:"), items, 0, false, &ok);
	if (ok && !item.isEmpty()) {
		for (int i=0; i<Plot::ncmax; i++) {
			if (graphResultSet[i] == 0) continue;
			if (item.compare(graphResultSet[i]->casename) == 0) {
				return i;
			}
		}
	}
	return -1;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::removeGraph()
{
	int i = selectGraphCase();
	if (i == -1) return;
	RESULT_SET *R = graphResultSet[i];
	// First remove the curves
	graph_nborn->removeCurve(R->casename);
	graph_ntot->removeCurve(R->casename);
	graph_nDC->removeCurve(R->casename);
	graph_teffgen->removeCurve(R->casename);
	graph_ncog->removeCurve(R->casename);
	if (!show_outputdata)
		graph_ncogseed->removeCurve(R->casename);
	// Then remove the graph case
	graphResultSet[i] = 0;
	nGraphCases--;
	drawGraphs();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::removeAllGraphs()
{
	clearAllGraphs();
}

//---------------------------------------------------------------------
// Updates input parameter (QLineEdit) widgets according to the slider
//---------------------------------------------------------------------
void MainWindow::updateSliderBox()
{
    paramSaved = false;          // Keeps track of the fact that a param has been changed but not saved.
    // ---Get value from slider
    QString slider_str = sender()->objectName();
    QString stag = slider_str.mid(7);
    int ival = ((QSlider *)sender())->value();

    // --- Get param index from workingParameterList
	int k;
	for (k=0; k<nParams; k++) {
		PARAM_SET p = parm->get_param(k);
		if (stag.compare(p.tag) == 0)
			break;
	}
    int j = param_to_sliderIndex[k];
    SliderPlus *sp = sliderplus_list[j];
    double v = sp->int_to_val(ival);
    QString vstr = sp->val_to_str(v);
    ((QLineEdit *)sliderParam[j])->setText(vstr);
}

//------------------------------------------------------------------------------------------------------
// changeParam() is invoked in response to a signal sent when a value in a QLineEdit etc. widget
// is changed.  Note that when a QRadioButton widget is changed, signals are sent both the radiobuttons
// that change, but only one signal is used to change the parameter value.
//------------------------------------------------------------------------------------------------------
void MainWindow::changeParam()
{
	paramSaved = false;
    QObject *w = sender(); // Gets the pointer to the object that invoked the changeParam slot.
	if (w->isWidgetType()) {
		QString wname = w->objectName();
		if (wname.contains("line_")) {
			QString wtag = wname.mid(5);
			QLineEdit *lineEdit = (QLineEdit *)w;
            QString text = lineEdit->displayText();
			// Determine if there is a slider associated with the sender widget
			for (int k=0; k<parm->nParams; k++) {
				PARAM_SET p = parm->get_param(k);
				if (wtag.compare(p.tag) == 0) {
					int j = param_to_sliderIndex[k];
					if (j >= 0) {
                        QSlider *slider = slider_list[j];
                        SliderPlus *sp = sliderplus_list[j];
                        double v = sp->str_to_val(text);
                        int ival = sp->val_to_int(v);
                        int ival_old = sp->val_to_int(p.value);
						if (ival != ival_old) {
                            slider->setSliderPosition(ival);
						}
					}
					parm->set_value(k,text.toDouble());
					break;
				}
			}
		} else if (wname.contains("spin_")) {
			QSpinBox *spinBox = (QSpinBox *)w;
            int v = spinBox->value();
			QString wtag = wname.mid(5);
			for (int k=0; k<parm->nParams; k++) {
				PARAM_SET p = parm->get_param(k);
				if (wtag.compare(p.tag) == 0) {
					parm->set_value(k,v);
                    break;
				}
				if (wname.contains("NCPU")) {
					ncpu = v;
				}
			}
		} else if (wname.contains("cbox_")) {

			QCheckBox *checkBox = (QCheckBox *)w;
			int v;
			QString wtag = wname.mid(5);
			for (int k=0; k<parm->nParams; k++) {
				PARAM_SET p = parm->get_param(k);
				if (wtag.compare(p.tag) == 0) {
					parm->set_value(k,v);
                    break;
				}
			}
			if (wname.contains("savepos")) {
				if (checkBox->isChecked()) {
					setSavePosStart();
				}
			}
		} else if (wname.contains("rbut_")) {
			QRadioButton *radioButton = (QRadioButton *)w;
			if (radioButton->isChecked()) {
				QString wtag = wname.mid(5);
				int rbutton_case;
				wtag = parse_rbutton(wtag,&rbutton_case);
				for (int k=0; k<parm->nParams; k++) {
					PARAM_SET p = parm->get_param(k);
					if (wtag.compare(p.tag) == 0) {
						parm->set_value(k,rbutton_case);
						break;
					}
				}
			}
		}
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::redrawDistPlot()
{
        int i_m = 0, i_s = 0;
    QString sname = sender()->objectName();
	for (int k=0; k<3; k++) {
		QwtPlot *qp = distplot_list[k];
        QString tag = qp->objectName().mid(8);
        QString tag_m = tag + "_MEDIAN";
        QString tag_s = tag + "_SHAPE";
		if (sname.endsWith(tag_m) || sname.endsWith(tag_s)) {   
			for (int i=0; i<nWidgets; i++) {
				QString wname = widget_list[i]->objectName();
                if (wname.endsWith(tag_m))
                    i_m = i;
                else if (wname.endsWith(tag_s))
                    i_s = i;
			}

            QString median_str = ((QLineEdit *)widget_list[i_m])->text();
            QString shape_str = ((QLineEdit *)widget_list[i_s])->text() ;
			if (median_str.compare("") == 0) return;
			if (shape_str.compare("") == 0) return;
            double median = median_str.toDouble();
            median = max(0.001, median);
            double shape = shape_str.toDouble();
            shape = max(1.001, shape);
			
			double *x = new double[nDistPts];
			double *prob = new double[nDistPts];
            create_lognorm_dist(median,shape,nDistPts,x,prob);
            int n = dist_limit(prob,nDistPts);
            double xmax = x[n];
            curve_list[k]->setData(x, prob, n);
            qp->setAxisScale(QwtPlot::xBottom, 0.0, xmax, 0.0);
            qp->replot();
			delete [] x;
			delete [] prob;
			x = 0;
			prob = 0;
		}
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
int MainWindow::dist_limit(double *p, int n)
{
	int i, imax;
    double pmax = 0;

	imax = n-1;
	for (i=0; i<n; i++) {
		if (p[i] > pmax) {
			pmax = p[i];
			imax = i;
		}
	}
    double plim = 0.01*pmax;
	for (i=n-1; i>0; i--) {
		if (p[i] > plim) {
            return min(n-1,max(i,2*imax));
		}
	}
	return 1;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
double MainWindow::erf(double z)
{
   double t = 1.0 / (1.0 + 0.5 * fabs(z));
   // use Horner's method
   double ans = 1 - t * exp( -z*z -  1.26551223 +
         t * ( 1.00002368 +
         t * ( 0.37409196 +
         t * ( 0.09678418 +
         t * (-0.18628806 +
         t * ( 0.27886807 +
         t * (-1.13520398 +
         t * ( 1.48851587 +
         t * (-0.82215223 +
         t * ( 0.17087277))))))))));
   if (z >= 0.0)
       return ans;
   else
       return -ans;
}

//-----------------------------------------------------------------------------------------
// When X is distributed N(mu,sig), this gives Prob{x1 < X <= x2}
//-----------------------------------------------------------------------------------------
double MainWindow::pnorm(double x1, double x2, double mu, double sig)
{
    double z1, z2, e1, e2;
		
	z1 = (x1-mu)/sig;
    e1 = erf(z1/sqrt(2.0))/2;
    z2 = (x2-mu)/sig;
    e2 = erf(z2/sqrt(2.0))/2;
    return e2 - e1;
}

//-----------------------------------------------------------------------------------------
// When log(X) is distributed N(mu,sig), this gives Prob{x1 < X <= x2}
//-----------------------------------------------------------------------------------------
double MainWindow::plognorm(double x1, double x2, double mu, double sig)
{
    double z1, z2, e1, e2;

    z1 = 0;
    z2 = 0;
    if (x1 == 0)
        e1 = -0.5;
	else {
        z1 = (log(x1)-mu)/sig;
        e1 = erf(z1/sqrt(2.0))/2;
	}
    if (x2 == 0)
        e2 = -0.5;
	else {
        z2 = (log(x2)-mu)/sig;
        e2 = erf(z2/sqrt(2.0))/2;
	}
    return e2 - e1;
}

//-----------------------------------------------------------------------------------------
// Create the lognormal distribution with median = p1, shape = p2
// at n points stored in x[], probability values stored in prob[].
// Note that x[0] = 0.
// The range of x is currently just less than 4*median.  This should be
// OK for values of shape < 2.
// Convert probability into probability density
//-----------------------------------------------------------------------------------------
void MainWindow::create_lognorm_dist(double p1, double p2,int n, double *x, double *prob)
{
	double xmax, dx, mu_l, sig_l, x1, x2;

    if (p1 >= 0.5)
        xmax = p1*4;
    else
        xmax = p1*8;
        
    dx = xmax/n;
    mu_l = log(p1);
    sig_l = log(p2);
	for (int ix=0; ix<n; ix++) {
        x1 = (ix - 0.5)*dx;
        x2 = x1 + dx;
		x1 = max(x1,0.0);
        x[ix] = (x1+x2)/2;
        prob[ix] = plognorm(x1,x2,mu_l,sig_l)/(x2-x1);
	}
}


//======================================================================================================
//------------------------------------------------------------------------------------------------------
// SliderPlus class member definitions
//------------------------------------------------------------------------------------------------------
SliderPlus::SliderPlus(QString aname, double valmin, double valmax, int nval, int iparam, int kwidget)
{
	int i;
    name = aname;
    pindex = iparam;
    windex = kwidget;
    dv = (valmax - valmin)/nval;
	for (i=10; i>-10; i--) {
		if (pow(10.0,i) < dv) {
            dv = pow(10.0,i);
            break;
		}
	}
    i = int(valmin/dv);
    vmin = dv*(i+1);
    int n1 = (int)((valmax - vmin)/dv + 0.5);	//round
	if (n1 > 5*nval) {
        dv = 5*dv;
        n1 = n1/5;
	}
	else if (n1 > 2*nval) {
        dv = 2*dv;
        n1 = n1/2;
	}
    i = int(valmin/dv);
    vmin = dv*i;
    n = n1;
    vmax = vmin + n*dv;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
SliderPlus::~SliderPlus()
{}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
int SliderPlus::val_to_int(double v) {
    int i = (int)((v-vmin)/dv + 0.5);	//round
    return i;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
double SliderPlus::int_to_val(int i) {
    double v = vmin + i*dv;
    return v;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
QString SliderPlus::val_to_str(double v) {
	int i = SliderPlus::val_to_int(v);
	QString vstr = QString::number(int_to_val(i));
    return vstr;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
double SliderPlus::str_to_val(QString vstr) {
    double v = vstr.toDouble();
    return v;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
int SliderPlus::pIndex() {
    return pindex;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
int SliderPlus::wIndex() {
    return windex;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
int SliderPlus::nTicks() {
    return n;
}
