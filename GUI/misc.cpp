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
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void SocketHandler::run()
{
//	QObject::moveToThread(this);
	sprintf(msg,"run: port: %d", port);
	LOG_MSG(msg);
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
	LOG_MSG("In processor");
    socket = tcpServer->nextPendingConnection();
	sprintf(msg,"got server connection: %p",socket);	
	LOG_MSG(msg);
    emit sh_connected();
	QString qdata;
	QByteArray ba;
	ba.resize(1024);
	while (true) {
		if (stopped) break;
//		if (port == CPORT0) {
//			LOG_MSG("CPORT0 waiting for ReadyRead");
//		} else {
//			LOG_MSG("CPORT1 waiting for ReadyRead");
//		}
		socket->waitForReadyRead(100);
//		if (port == CPORT0) {
//			LOG_MSG("CPORT0 is ReadyRead");
//		} else {
//			LOG_MSG("CPORT1 is ReadyRead");
//		}
		int nb = socket->bytesAvailable();
		if (nb > 0) {
			ba = socket->readLine(1024);
			qdata = QString(ba);
//		    emit sh_output(qdata); // Emit signal to update GUI
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
			} else {
//				LOG_MSG("No bytes yet");
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
ExecThread::ExecThread(QString infile, QString dllfile)
{
	inputFile = infile;
	dll_path = dllfile;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::run()
{
	LOG_MSG("Invoking DLL...");
	const char *infile, *outfile, *resfile, *runfile;
	QString infile_path = inputFile;
	QString casename = QFileInfo(inputFile).baseName();
    QString outfile_path = casename + ".res";	
    QString resfile_path = "res.log";
    QString runfile_path = "run.log";

		LOG_QMSG(infile_path);

	int len_infile = infile_path.length();
	int len_outfile = outfile_path.length();
	int len_resfile = resfile_path.length();
	int len_runfile = runfile_path.length();
	std::string std_infile = infile_path.toStdString();
	std::string std_outfile = outfile_path.toStdString();
	std::string std_resfile = resfile_path.toStdString();
	std::string std_runfile = runfile_path.toStdString();
	infile = std_infile.c_str();
	outfile = std_outfile.c_str();
	resfile = std_resfile.c_str();
	runfile = std_runfile.c_str();

#ifdef __COMPILETIME_LOADING__

#ifdef __GFORTRAN_DLL__
	__bone_mod_MOD_execute(&nsteps);
//					const_cast<char *>(infile),   \
//					const_cast<char *>(outfile),  \
//					const_cast<char *>(resfile),  \
//					const_cast<char *>(runfile),  \
//					len_infile,  \
//					len_outfile, \
//					len_resfile, \
//					len_runfile);
#else
	paused = false;
	execute(const_cast<char *>(infile),&len_infile);
	LOG_MSG("execute returned");
	get_dimensions(&NX,&NY,&NZ,&NBY);
	sprintf(msg,"exthread: nsteps: %d",nsteps);
	LOG_MSG(msg);
	for (int i=1; i<= nsteps; i++) {
		bool updated = false;
		if (paused && !updated) {
			snapshot();
			updated = true;
		}
		while(paused) {
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
	}
	snapshot();
	Sleep(10);
	int result = 0;
	terminate_run(&result);
#endif
	LOG_MSG("Returning from exthread");
	return;

#else

//	sprintf(msg,"inputFile: %p  %p",(const char *)((infile_path.toStdString()).data()), infile);
//	LOG_MSG(msg);
	LOG_QMSG(outfile_path);

	QLibrary myLib(dll_path);

	LOG_QMSG(dll_path);

	LOG_QMSG(myLib.errorString());

#ifdef __GFORTRAN_DLL__
//	typedef int (*MyPrototype)(int *, char *, char *, char *, char *, int, int, int, int);
//	typedef int (*MyPrototype)(int *);
	typedef int (*MyPrototype)(char *, int);
	LOG_MSG("__GNUC__ is defined");
	MyPrototype execute = (MyPrototype) myLib.resolve("__bone_mod_MOD_execute");	// NOTE: DLL procedure name must be in uppercase
	if (execute) {
			LOG_MSG("resolved EXECUTE()");
			execute(const_cast<char *>(infile), len_infile);
//			execute(&nsteps);
//					const_cast<char *>(infile),   \
//					const_cast<char *>(outfile),  \
//					const_cast<char *>(resfile),  \
//					const_cast<char *>(runfile),  \
//					len_infile,  \
//					len_outfile, \
//					len_resfile, \
//					len_runfile);

#else

//	typedef int (*MyPrototype)(int *, char *, int, char *, int, char *, int, char *, int);
//	typedef int (*MyPrototype)(int *);
	typedef int (*MyPrototype)(char *, int);
	MyPrototype execute = (MyPrototype) myLib.resolve("EXECUTE");	// NOTE: DLL procedure name must be in uppercase
	if (execute) {	
		LOG_MSG("resolved EXECUTE()");
		execute(const_cast<char *>(infile), len_infile);
//		execute(&nsteps);
//			const_cast<char *>(infile), len_infile, \
//			const_cast<char *>(outfile), len_outfile, \
//			const_cast<char *>(resfile), len_resfile, \
//			const_cast<char *>(runfile), len_runfile);
#endif

	    LOG_MSG("DLL execution completed");
//		sleep(1);
//		LOG_MSG("unload DLL");
		myLib.unload();
//		LOG_MSG("deleteLater");
		myLib.deleteLater();
//		LOG_MSG("end run");
			} else {
				LOG_MSG("Failed to resolve EXECUTE in the dynamic library");
				LOG_QMSG(myLib.errorString());
				exit(1);
	}

//	done = true;
#endif
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::snapshot()
{
//	LOG_MSG("snapshot");
	mutex2.lock();
	get_scene(&ncap_list,cap_list,&nmono_list,mono_list,&npit_list,pit_list,&nclast_list,clast_list);
//	sprintf(msg,"got: ncap %d nmono %d npit %d nclast %d",ncap_list,nmono_list,npit_list,nclast_list);
//	LOG_MSG(msg);
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
