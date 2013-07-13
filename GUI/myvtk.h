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
//#include <vtkMPEG2Writer.h>
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
#include <vtkTextSource.h>
#include <vtkGlyph3D.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPolyDataMapper.h>
#include <vtkUnsignedCharArray.h>
#include <vtkDoubleArray.h>


using namespace std;

struct cell_pos {
	int tag;
	int x, y, z;
	int state;
};
typedef cell_pos CELL_POS;

struct bond_pos {
	int TCtag;
	int DCtag;
};
typedef bond_pos BOND_POS;

struct capillary_segment {
	double pos1[3], pos2[3];
	double radius;
};
typedef capillary_segment CAPILLARY_SEGMENT;

struct pit_pos {
	int pos[3];
	double ypit;
};
typedef pit_pos PIT_POS;

struct surface_tile {
	int pos[3];
	char axis;
	bool pit;
};
typedef surface_tile SURFACE_TILE;

struct bone {
	double y;
	vtkActor *actor;
};
typedef bone BONE;

struct vert_tile {
	double pos[3];
	double height;
	char axis;
	vtkActor *actor;
};
typedef vert_tile VERT_TILE;

struct clast_pos {
	double pos[3];
	double dir[3];
	double size;
};
typedef clast_pos CLAST_POS;

struct blast_pos {
	double pos[3];
	int state;
};
typedef blast_pos BLAST_POS;

class MyVTK
{
public:
	MyVTK(QWidget *);
	~MyVTK();

	void read_cell_positions(QString, QString, bool);
	void get_cell_positions(bool fast);
	void init();
	void cleanup();
	void renderCells();
	void makeBox();
	void process_Mcells();
	void process_Tcells();
    void process_Dcells(bool);
    void process_bonds();
	void process_capillaries();
	void oldprocess_tiles();
	void process_tiles();
	void oldcreateTileList();
	void createTileList();
	void process_osteoclasts();
	void process_osteoblasts();
	int inTileList(SURFACE_TILE *tile);
	bool startPlayer(QString, QTimer *, bool);
	bool nextFrame();
	void pause();
	void playon();
	void saveSnapshot(QString, QString);
	void stop();
	void makeEllipsoid();

	QList<CELL_POS > MCpos_list;
	QList<CELL_POS > TCpos_list;
	QList<CELL_POS > DCpos_list;
	QList<BOND_POS > bondpos_list;
	QList<CAPILLARY_SEGMENT> capillary_list;
	QList<PIT_POS > pitpos_list;
	QList<CLAST_POS > clastpos_list;
	QList<BLAST_POS > blastpos_list;
	QList<SURFACE_TILE> tile_list;
	QList<VERT_TILE> verttile_list;
	QList<vtkActor *> M_Actor_list;
	QList<vtkActor *> T_Actor_list;
	QList<vtkActor *> D_Actor_list;
	QList<vtkActor *> B_Actor_list;
	QList<vtkActor *> C_Actor_list;
	QList<vtkActor *> OC_Actor_list;
	QList<vtkActor *> OB_Actor_list;
	QList<vtkActor *> Tile_Actor_list;

	BONE **bone_array;

	QVTKWidget* qvtkWidget;
	vtkRenderWindow *renWin;	
	vtkRenderer* ren;
	vtkRenderWindowInteractor * iren;
	vtkPolyDataMapper *McellMapper;
	vtkPolyDataMapper *TcellMapper;
	vtkPolyDataMapper *DcellMapper;
	vtkPolyDataMapper *bondMapper;
	vtkPolyDataMapper *capillaryMapper;
	vtkPolyDataMapper *cylMapper;
	vtkPolyDataMapper *tileMapper;
	vtkPolyDataMapper *clastMapper;
	vtkPolyDataMapper *blastMapper;
	vtkPolyDataMapper *textMapper;

	vtkSmartPointer<vtkPolyData> polygonPolyData;
	vtkSmartPointer<vtkPoints> points;
	vtkSmartPointer<vtkUnsignedCharArray> colors;
	vtkSmartPointer<vtkPolyData> ginput;
	vtkPolyDataMapper *squareMapper;
	vtkSmartPointer<vtkGlyph3D> glypher;
	vtkActor *glactor;

//	vtkMPEG2Writer *mpg;
//	vtkSmartPointer<vtkPNGWriter> writer;
//	vtkSmartPointer<vtkBMPWriter> writer;
	vtkSmartPointer<vtkJPEGWriter> writer;
//	vtkSmartPointer<vtkTIFFWriter> writer;
	vtkSmartPointer<vtkImageCast> castFilter;
	vtkWindowToImageFilter *w2i;
	vtkWindowToImageFilter *w2img;

	char msg[2048];
	double zoomlevel;
	double Pi;
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
	vtkActor *boxActor[8];

};

#endif
