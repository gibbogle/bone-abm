#include <string>
#include <fstream>
#ifdef _WIN32
#include <windows.h>
#endif
#include <QTcpServer>
#include <QTcpSocket>
#include <QtGui>
#include <QTcpServer>
#include <QMessageBox>

#include "misc.h"
#include "log.h"
#include "transfer.h"

#include "libbone.h"

LOG_USE();
char msg[2048];

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
SocketHandler::SocketHandler(int newport, QObject *parent)
	: QThread(parent)
{
    exiting = false;
	stopped = false;
    port = newport;
	sprintf(msg,"SocketHandler: port: %d",port);
	LOG_MSG(msg);
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
SocketHandler::~SocketHandler() // make sure the worker object is destroyed   
{
    exiting = true;
    wait();
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void SocketHandler::stop()
{
	stopped = true;
    LOG_MSG("stopped in stop");
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void SocketHandler::run()
{
	quint16 qport = port;
	QString addressStr = "127.0.0.1";
	QHostAddress hostAddress;
	hostAddress.setAddress(addressStr);
    tcpServer = new QTcpServer(this);
	connect(tcpServer, SIGNAL(newConnection()), this, SLOT(processor()), Qt::DirectConnection);
    if (!tcpServer->listen(hostAddress,qport)) {
		sprintf(msg,"Unable to start the server: port: %d", port);
		LOG_MSG(msg);
        return;
    }
	sprintf(msg,"Listening on port: %d",tcpServer->serverPort());
	LOG_MSG(msg);
	LOG_MSG("serverAddress:");
	LOG_QMSG((tcpServer->serverAddress()).toString());
	bool timedOut = false;
	tcpServer->waitForNewConnection(-1,&timedOut);
	exec();
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void SocketHandler::processor()
{
    socket = tcpServer->nextPendingConnection();
	sprintf(msg,"got server connection: %p",socket);	
	LOG_MSG(msg);
    emit sh_connected();
	QString qdata;
	QByteArray ba;
	ba.resize(1024);
	while (true) {
		if (stopped) break;
		socket->waitForReadyRead(100);
		int nb = socket->bytesAvailable();
		if (nb > 0) {
			ba = socket->readLine(1024);
			qdata = QString(ba);
			QStringList s = qdata.split("^",QString::SkipEmptyParts);
			for (int k=0; k<s.length(); k++) {
				emit sh_output(s[k]); // Emit signal to update GUI
				if (port == CPORT0) {
					LOG_QMSG(s[k]);
				}
			}
			if (quitMessage(qdata)) {
				sprintf(msg,"Closing connection: port: %d", port);
				LOG_MSG(msg);
		        break;
			}
		}
	}
	socket->close();
	tcpServer->close();
	if (port == CPORT0) {
		emit sh_disconnected();		// Is it right that both threads emit this?
		LOG_MSG("emitted sh_disconnected");
	}
}


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
ExecThread::ExecThread(int runcase, QString infile)
{
    run_case = runcase;
	inputFile = infile;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::run()
{
    int result;
	const char *infile;
	QString infile_path = inputFile;
	QString casename = QFileInfo(inputFile).baseName();
	int len_infile = infile_path.length();
	std::string std_infile = infile_path.toStdString();
	infile = std_infile.c_str();

	paused = false;
    execute(const_cast<char *>(infile),&len_infile, &run_case, &result);
    if (result != 0) {
        terminate_run(&result);
        return;
    }
	// Now the "Stop" button can be enabled
	emit(initialized());
	get_dimensions(&NX,&NY,&NZ,&NBY);
	sprintf(msg,"exthread: nsteps: %d",nsteps);
	LOG_MSG(msg);
	for (int i=1; i<= nsteps; i++) {
		bool updated = false;
		if (paused && !updated) {
			snapshot();
			updated = true;
		}
		while(paused || leftb) {
			Sleep(100);
		}
		if (stopped) break;
		int res=0;
		simulate_step(&res);
		if (stopped) break;
		if (i%240 == 0) {
			mutex1.lock();
			get_summary(summaryData);
			mutex1.unlock();
			emit summary();		// Emit signal to update summary plots, at hourly intervals
		}
		if (stopped) break;
		if (i%nt_vtk == 0) {
			if (showingVTK != 0) {
				snapshot();
				Sleep(10);
			}
		}
		if (stopped) break;
		if (res == 1) break;
	}
	snapshot();
	Sleep(10);
	terminate_run(&result);
	LOG_MSG("Returning from exthread");
	return;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::snapshot()
{
	mutex2.lock();
	get_scene(&ncap_list,cap_list,&nmono_list,mono_list,&npit_list,pit_list,&nclast_list,clast_list,&nblast_list,blast_list);
	if (ncap_list > MAX_CAP) {
		LOG_MSG("Error: MAX_CAP exceeded");
		exit(1);
	}
	if (nmono_list > MAX_MONO) {
		LOG_MSG("Error: MAX_MONO exceeded");
		exit(1);
	}
	if (npit_list > MAX_PIT) {
		LOG_MSG("Error: MAX_PIT exceeded");
		exit(1);
	}
	if (nclast_list > MAX_CLAST) {
		LOG_MSG("Error: MAX_CLAST exceeded");
		exit(1);
	}
	if (nblast_list > MAX_BLAST) {
		LOG_MSG("Error: MAX_BLAST exceeded");
		exit(1);
	}
	mutex2.unlock();
	emit display(); // Emit signal to update VTK display
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::stop()
{
	stopped = true;
}

	//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::pause()
{
	paused = true;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::unpause()
{
	paused = false;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool quitMessage(QString msg)
{
	if (msg.contains("__EXIT__",Qt::CaseSensitive))
		return true;
	else
		return false;
}
