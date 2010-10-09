// myvtk.h
#ifndef MYVTK_H
#define MYVTK_H

#include <QtGui>
#include <QtCore>
#include <QIODevice>
#include <QVTKWidget.h>
#include <vtkRenderer.h> 
#include <vtkRenderWindow.h>
#include "vtkSphereSource.h"
#include "vtkCylinderSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include <vtkMPEG2Writer.h>
#include <vtkPNGWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkTIFFWriter.h>
#include <vtkBMPWriter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkSmartPointer.h>
#include <vtkImageCast.h>
#include <vtkPolygon.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>

//#include <vtkConfigure.h>

using namespace std;

struct cell_pos {
	int tag;
	int x, y, z;
	double diameter;
	double state;
};
typedef cell_pos CELL_POS;

struct bond_pos {
	int TCtag;
	int DCtag;
};
typedef bond_pos BOND_POS;

struct capillary_segment {
	int pos1[3], pos2[3];
	double radius;
};
typedef capillary_segment CAPILLARY_SEGMENT;

struct pit_pos {
	int pos[3];
};
typedef pit_pos PIT_POS;

struct surface_tile {
	int pos[3];
	char axis;
};
typedef surface_tile SURFACE_TILE;

class MyVTK
{
public:
	MyVTK(QWidget *);
	~MyVTK();

	void read_cell_positions(QString, QString, bool);
	void init();
	void cleanup();
	void renderCells(bool,bool);
	void process_Mcells();
	void process_Tcells();
    void process_Dcells(bool);
    void process_bonds();
	void process_capillaries();
	void process_tiles();
	void createTileList();
	int inTileList(SURFACE_TILE *tile);
	bool startPlayer(QString, QTimer *, bool);
	bool nextFrame();
	void pause();
	void playon();
	void saveSnapshot(QString, QString);
	void stop();

	QList<CELL_POS > MCpos_list;
	QList<CELL_POS > TCpos_list;
	QList<CELL_POS > DCpos_list;
	QList<BOND_POS > bondpos_list;
	QList<CAPILLARY_SEGMENT> capillary_list;
	QList<PIT_POS > pitpos_list;
	QList<SURFACE_TILE> tile_list;
	QList<vtkActor *> M_Actor_list;
	QList<vtkActor *> T_Actor_list;
	QList<vtkActor *> D_Actor_list;
	QList<vtkActor *> B_Actor_list;
	QList<vtkActor *> C_Actor_list;
	QList<vtkActor *> Tile_Actor_list;

	QVTKWidget* qvtkWidget;
	vtkRenderWindow *renWin;	
	vtkRenderer* ren;
	vtkRenderWindowInteractor * iren;
	vtkPolyDataMapper *McellMapper;
	vtkPolyDataMapper *TcellMapper;
	vtkPolyDataMapper *DcellMapper;
	vtkPolyDataMapper *bondMapper;
	vtkPolyDataMapper *capillaryMapper;
	vtkPolyDataMapper *tileMapper;

	vtkMPEG2Writer *mpg;
//	vtkSmartPointer<vtkPNGWriter> writer;
//	vtkSmartPointer<vtkBMPWriter> writer;
	vtkSmartPointer<vtkJPEGWriter> writer;
//	vtkSmartPointer<vtkTIFFWriter> writer;
	vtkSmartPointer<vtkImageCast> castFilter;
	vtkWindowToImageFilter *w2i;
	vtkWindowToImageFilter *w2img;
//	vtkSmartPointer<vtkPNGWriter> pngwriter;
//	vtkSmartPointer<vtkJPEGWriter> jpgwriter;

	char msg[2048];
	double zoomlevel;
	bool DCmotion;
	bool DCfade;
	bool first_VTK;
	bool playing;
	bool paused;
	bool save_image;
	QString casename;
	int framenum;
	QTimer *timer;
	QString infile;
	QFile *playerData;
	QTextStream *playerStream;

//	double *TCColor, *DCColor;

//	vtkSmartPointer<vtkPoints> points;
//	vtkSmartPointer<vtkCellArray> vertices;
//	vtkSmartPointer<vtkPolygon> polygon;
//	vtkSmartPointer<vtkCellArray> polygons;
//	vtkSmartPointer<vtkPolyData> polygonPolyData;
};

#endif