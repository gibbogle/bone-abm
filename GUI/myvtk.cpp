// myvtk.cpp

#include <vtkCamera.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>


#ifdef _WIN32
#include "windows.h"
#endif
#include "myvtk.h"
#include "log.h"
#include "transfer.h"

LOG_USE();


// Define interaction style
class MouseInteractorStyle4 : public vtkInteractorStyleTrackballCamera
{
  public:
	static MouseInteractorStyle4* New();
	vtkTypeMacro(MouseInteractorStyle4, vtkInteractorStyleTrackballCamera);

	virtual void OnLeftButtonDown()
	{
	  std::cout << "Pressed left mouse button." << std::endl;
	  LOG_QMSG("Pressed left mouse button.");
	  leftb = true;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
	}

	virtual void OnMiddleButtonDown()
	{
	  std::cout << "Pressed middle mouse button." << std::endl;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnMiddleButtonDown();
	}

	virtual void OnRightButtonDown()
	{
	  std::cout << "Pressed right mouse button." << std::endl;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnRightButtonDown();
	}

	virtual void OnLeftButtonUp()
	{
	  std::cout << "Released left mouse button." << std::endl;
	  LOG_QMSG("Released left mouse button.");
	  leftb = false;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnLeftButtonUp();
	}

	virtual void OnMiddleButtonUp()
	{
	  std::cout << "Released middle mouse button." << std::endl;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnMiddleButtonUp();
	}

	virtual void OnRightButtonUp()
	{
	  std::cout << "Released right mouse button." << std::endl;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnRightButtonUp();
	}

};

vtkStandardNewMacro(MouseInteractorStyle4);


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
MyVTK::MyVTK(QWidget *page)
{
	zoomlevel = 0.7;
	double backgroundColor[] = {0.0,0.0,0.0};

	Pi = 4*atan(1.0);
	leftb = false;
	qvtkWidget = new QVTKWidget(page,QFlag(0));
	LOG_MSG("Created a new QVTKWidget");
	QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(qvtkWidget);

	// Associate the layout with page_VTK
    page->setLayout(layout);

	// Create a renderer, and add it to qvtkWidget's render window.
	// The renderer renders into the render window. 
	ren = vtkRenderer::New();     
    renWin = qvtkWidget->GetRenderWindow();
    renWin->AddRenderer(ren);
	ren->SetBackground(backgroundColor);
//	ren->SetBackground(0.1, 0.2, 0.4);		// backgroundColor
	ren->ResetCamera();
	iren = qvtkWidget->GetInteractor();
//	iren->Initialize();
//	ren->RemoveAllViewProps();

	vtkSmartPointer<MouseInteractorStyle4> style =
	  vtkSmartPointer<MouseInteractorStyle4>::New();
	iren->SetInteractorStyle( style );

	iren->Initialize();

	// Create mappers
	vtkSphereSource *Tcell = vtkSphereSource::New();
    Tcell->SetThetaResolution(12);
    Tcell->SetPhiResolution(12);
    Tcell->SetRadius(0.5);
	TcellMapper = vtkPolyDataMapper::New();
	TcellMapper->SetInputConnection(Tcell->GetOutputPort());

	vtkSphereSource *Dcell = vtkSphereSource::New();
    Dcell->SetThetaResolution(12);
    Dcell->SetPhiResolution(12);
    Dcell->SetRadius(1.0);
	DcellMapper = vtkPolyDataMapper::New();
	DcellMapper->SetInputConnection(Dcell->GetOutputPort());

	vtkSphereSource *Mcell = vtkSphereSource::New();
	Mcell->SetThetaResolution(12);
	Mcell->SetPhiResolution(12);
	Mcell->SetRadius(0.5);
	McellMapper = vtkPolyDataMapper::New();
	McellMapper->SetInputConnection(Mcell->GetOutputPort());

	vtkCylinderSource *cyl = vtkCylinderSource::New();
	cyl->SetResolution(12);
	cyl->SetRadius(0.2);
	cyl->SetHeight(1);
	cylMapper = vtkPolyDataMapper::New();
	cylMapper->SetInputConnection(cyl->GetOutputPort());

	vtkCylinderSource *capillary = vtkCylinderSource::New();
	capillary->SetResolution(12);
	capillary->SetRadius(1);
	capillary->SetHeight(1);
	capillaryMapper = vtkPolyDataMapper::New();
	capillaryMapper->SetInputConnection(capillary->GetOutputPort());

	//Create the square
	// Setup points
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
//	vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
	points->InsertNextPoint(-0.5, 0.0, -0.5);
	points->InsertNextPoint( 0.5, 0.0, -0.5);
	points->InsertNextPoint( 0.5, 0.0,  0.5);
	points->InsertNextPoint(-0.5, 0.0,  0.5);
	// a unit square in the XY plane, at y=0

	/*
	vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
	polygon->GetPointIds()->SetNumberOfIds(4); //make a quad
	polygon->GetPointIds()->SetId(0, 0);
	polygon->GetPointIds()->SetId(1, 1);
	polygon->GetPointIds()->SetId(2, 2);
	polygon->GetPointIds()->SetId(3, 3);

	//Add the polygon to a list of polygons
	vtkSmartPointer<vtkCellArray> polygons = vtkSmartPointer<vtkCellArray>::New();
	polygons->InsertNextCell(polygon);

	//Create a PolyData
	vtkSmartPointer<vtkPolyData> polygonPolyData = vtkSmartPointer<vtkPolyData>::New();
	polygonPolyData->SetPoints(points);
	polygonPolyData->SetPolys(polygons);
	*/

	vtkSmartPointer<vtkIdList> pointsIds= vtkSmartPointer< vtkIdList>::New();

	for (int i=0;i<4;i++) pointsIds->InsertNextId(i);

	//Create a PolyData

	polygonPolyData = vtkSmartPointer<vtkPolyData>::New();

	polygonPolyData->Allocate();

	polygonPolyData->SetPoints(points);

	polygonPolyData->InsertNextCell(VTK_QUAD, pointsIds);


	//Create a mapper
	tileMapper = vtkPolyDataMapper::New();
	tileMapper->SetInput(polygonPolyData);
	tileMapper->ScalarVisibilityOff();

	vtkSmartPointer<vtkTextSource> textSource = vtkSmartPointer<vtkTextSource>::New();
	textSource->SetText("Hello");
	textSource->SetForegroundColor(1.0, 1.0, 1.0);
	textSource->BackingOff();
	textSource->Update();

	// Create a mapper
	textMapper = vtkPolyDataMapper::New();
	textMapper->SetInputConnection(textSource->GetOutputPort());

	// Create a glypher
//	Even though a programmable source like that wiki example is probably the "cleanest" way
//	to do things, there is a simpler way to build a dataset for vtkGlyph3D.

//	The vtkGlyph3D filter needs two inputs: one input for the positions of each glyph
//	(i.e. one point per glyph) and another input for the glyph shape
//	(i.e. the square that you already made).

//	The first input is just a list of points, you can build it as a polydata just like
//	you did with the squares, except that you want to use "verts" as your cells:
	/*
	NX = 100; NZ = 100;
	int npos = NX*NZ;
	vtkSmartPointer<vtkPoints> positions = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> dummycells = vtkSmartPointer<vtkCellArray>::New();
	positions->SetNumberOfPoints(npos);         // number of squares
	dummycells->InsertNextCell(npos);           // size of dummy "verts" array
	int x, z, ipos=0;
	double pos[3];
	for (x = 1; x < NX+1; x++) {
		for (z = 1; z < NZ+1; z++) {
			pos[0] = x;
			pos[1] = NBY + 0.5;
			pos[2] = z;
			positions->SetPoint(ipos, pos[0], pos[1], pos[2]);     // square positions set here
			dummycells->InsertCellPoint(ipos);
			ipos++;
		}
	}
	vtkSmartPointer<vtkPolyData> input = vtkSmartPointer<vtkPolyData>::New();
	input->SetPoints(positions);
	input->SetVerts(dummycells);

	glypher = vtkSmartPointer<vtkGlyph3D>::New();
	glypher->SetInput(input);
	glypher->SetSource(polygonPolyData);
	glypher->SetScaleFactor(1.0);
	*/
//	Then feed the output of the glypher into the mapper, and it will automatically
//	draw a square at each position.
//	squareMapper = vtkPolyDataMapper::New();
	/*
	squareMapper->SetInputConnection(glypher->GetOutputPort());
	vtkActor *glactor = vtkActor::New();
	glactor->SetMapper(squareMapper);
	double boneColor[] = {0.9,0.9,0.5};
	glactor->GetProperty()->SetColor(boneColor);
	glactor->GetProperty()->SetAmbient(0.5);
	glactor->GetProperty()->SetDiffuse(0.2);
	glactor->GetProperty()->SetSpecular(0.5);
	ren->AddActor(glactor);
	*/
//	To change the positions, you would do this:
//	positions->SetPoint(idx, x, y, z); // set point idx to (x,y,z)
//	glypher->Modified(); // tell glypher that the data has changed

//	Then re-render, and the squares will be at their new positions. I have not tested the
//	code above, but I've done similar things in the past.
//  David Gobbi

	makeEllipsoid();

	//legendScaleActor = vtkSmartPointer<vtkLegendScaleActor>::New();

	// Create image filter for save Snapshot()
	w2img = vtkWindowToImageFilter::New();
//	pngwriter = vtkSmartPointer<vtkPNGWriter>::New();
//	jpgwriter = vtkSmartPointer<vtkJPEGWriter>::New();

	bone_array = 0;
	glypher = 0;
	squareMapper = 0;
	glactor = 0;
	first_VTK = true;
	DCmotion = false;
	DCfade = true;
	playing = false;
	paused = false;

	ren->GetActiveCamera()->Zoom(zoomlevel);		// try zooming OUT
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
MyVTK::~MyVTK()
{
}

//-----------------------------------------------------------------------------------------
// The ellipsoidal surface of the osteoclast is defined by three parameters.
//	a = length in the x direction
//	b = height in the y direction
//	c = width in the z direction
// The equation for the surface is:
//	x^2/a^2 + y^2/b^2 + z^2/c^2 = 1
// To discretize:
// Each quadrant is treated in the same way.  The projection of the quadrant onto the x-z
// plane is sliced into N segments (equal angles in the x-z projection), then the edges
// of the segments are divided into M pieces.  This generates, for each segment, M-1
// quads and 1 triangle.  The corresponding points on the surface are found from the eqtn.
//-----------------------------------------------------------------------------------------
void MyVTK::makeEllipsoid()
{
	static const int N = 5;	// number of segments
	static const int M = 6;	// number of divisions of a segment
	double a = 1.0;
	double b = 0.2;
	double c = 1.0;
	int pt_id[4*N+1][M+1];
	double beta[M+1];	// this is the fractional distance along the slice to a node
	int i, islice, id;

	// start with equal spacing
	for (i=0; i<M+1; i++) {
		beta[i] = i/double(M);
	}

	double dalpha = Pi/(2*N);
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	id = 0;
	for (islice = 0; islice<4*N+1; islice++) {
		double x, y, z, xm, zm, r, alpha, tana;
		alpha = islice*dalpha;
		// The line is given by z/x = tan(alpha), or z = x.tan(alpha)
		if (cos(alpha) == 0) {
			xm = 0;
			zm = c;
		} else {
			tana = sin(alpha)/cos(alpha);
			xm = sqrt(1./(1./(a*a) + tana*tana/(c*c)));
			zm = xm*tana;
		}
		r = sqrt(xm*xm+zm*zm);
		xm = r*cos(alpha);
		zm = r*sin(alpha);
		for (i=0; i<M+1; i++) {
			x = beta[i]*xm;
			z = beta[i]*zm;
			if (i == M)
				y = 0;
			else
				y = b*sqrt(1 - x*x/(a*a) - z*z/(c*c));
			if (id > 0 && i == 0) {
				pt_id[islice][i] = 0;
			} else {
				pt_id[islice][i] = id;
				id++;
				points->InsertNextPoint(x,y,z);	// # id
			}
		}
	}

	vtkSmartPointer<vtkCellArray> polygons = vtkSmartPointer<vtkCellArray>::New();
	for (islice = 0; islice<4*N; islice++) {
		// make a triangle
		vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
		polygon->GetPointIds()->SetNumberOfIds(3);
		polygon->GetPointIds()->SetId(0, 0);
		polygon->GetPointIds()->SetId(1, pt_id[islice][1]);
		polygon->GetPointIds()->SetId(2, pt_id[islice+1][1]);
		polygons->InsertNextCell(polygon);
		for (i=1; i<M; i++) {
			// make a quad
			vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
			polygon->GetPointIds()->SetNumberOfIds(4);
			polygon->GetPointIds()->SetId(0, pt_id[islice][i]);
			polygon->GetPointIds()->SetId(1, pt_id[islice][i+1]);
			polygon->GetPointIds()->SetId(2, pt_id[islice+1][i+1]);
			polygon->GetPointIds()->SetId(3, pt_id[islice+1][i]);
			polygons->InsertNextCell(polygon);
		}
	}

	//Create a PolyData
	vtkSmartPointer<vtkPolyData> polygonPolyData = vtkSmartPointer<vtkPolyData>::New();
	polygonPolyData->SetPoints(points);
	polygonPolyData->SetPolys(polygons);

	//Create a mapper
	clastMapper = vtkPolyDataMapper::New();
	clastMapper->SetInput(polygonPolyData);
	clastMapper->ScalarVisibilityOff();

}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::get_cell_positions(bool fast)
{
//	LOG_MSG("VTK: get_cell_positions");
	MCpos_list.clear();
	capillary_list.clear();
	pitpos_list.clear();
	clastpos_list.clear();
	int maxtag = 0;
	for (int i=0; i<nmono_list; i++) {
		int j = 5*i;
		CELL_POS cp;
		int status = mono_list[j+4];
		if (fast && status < 2) continue;
		cp.tag = mono_list[j];
		cp.x = mono_list[j+1];
		cp.y = mono_list[j+2];
		cp.z = mono_list[j+3];
		cp.state = status;
		MCpos_list.append(cp);
		if (cp.tag > maxtag) maxtag = cp.tag;
	}
//	sprintf(msg,"npit_list: %d",npit_list);
//	LOG_MSG(msg);
	for (int i=0; i<npit_list; i++) {
		int j = 4*i;
		PIT_POS pp;
		pp.pos[0] = pit_list[j];
		pp.pos[1] = pit_list[j+1];
		pp.pos[2] = pit_list[j+2];
		pp.ypit = pit_list[j+3];
		pitpos_list.append(pp);
	}
	for (int i=0; i<ncap_list; i++) {
		int j = 7*i;
		CAPILLARY_SEGMENT cs;
		cs.pos1[0] = cap_list[j];
		cs.pos1[1] = cap_list[j+1];
		cs.pos1[2] = cap_list[j+2];
		cs.pos2[0] = cap_list[j+3];
		cs.pos2[1] = cap_list[j+4];
		cs.pos2[2] = cap_list[j+5];
		cs.radius  = cap_list[j+6];
		capillary_list.append(cs);
	}
	for (int i=0; i<nclast_list; i++) {
		int j = 7*i;
		CLAST_POS oc;
		oc.pos[0] = clast_list[j];
		oc.pos[1] = clast_list[j+1];
		oc.pos[2] = clast_list[j+2];
		oc.dir[0] = clast_list[j+3];
		oc.dir[1] = clast_list[j+4];
		oc.dir[2] = clast_list[j+5];
		oc.size = clast_list[j+6];
		clastpos_list.append(oc);
	}
//	LOG_MSG("VTK: did get_cell_positions");
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::read_cell_positions(QString infileName, QString outfileName, bool savepos)
{	
//	LOG_MSG("VTK: read_cell_positions");
    TCpos_list.clear();
	MCpos_list.clear();
//	bondpos_list.clear();
	capillary_list.clear();
	pitpos_list.clear();
	QString line, saveline;
	QTextStream *out = NULL;
	QFile *vtkdata = NULL;
	if (savepos) {
		vtkdata = new QFile(outfileName);
		if (!vtkdata->open(QFile::Append )) {
			LOG_MSG("Open failure on vtk file");
			return;
		}
		out = new QTextStream(vtkdata);
	}
//	LOG_MSG("VTK: ready to read posdata");
	QFile posdata(infileName);
	if (posdata.open(QFile::ReadOnly)) {
		QTextStream in(&posdata);
		do {
			line = in.readLine();
			if (line.length() > 0) {
				if (savepos) {
					*out << line << "\n";
					out->flush();
				}
				saveline = line;
				QStringList s = line.split(" ",QString::SkipEmptyParts);
				if (s[0].compare("T") == 0) {
					CELL_POS cp;
					cp.tag = s[1].toInt();
					cp.x = s[2].toInt();
					cp.y = s[3].toInt();
					cp.z = s[4].toInt();
//					cp.diameter = s[5].toDouble();
//					cp.state = s[6].toDouble();
					TCpos_list.append(cp);
				} else if (s[0].compare("M") == 0) {
					CELL_POS cp;
					cp.tag = s[1].toInt();
					cp.x = s[2].toInt();
					cp.y = s[3].toInt();
					cp.z = s[4].toInt();
//					cp.diameter = s[5].toDouble();
					cp.state = s[5].toInt();
					MCpos_list.append(cp);
				} else if (s[0].compare("B") == 0) {
					NX = s[1].toInt();
					NY = s[2].toInt();
					NZ = s[3].toInt();
					NBY = s[4].toInt();
				} else if (s[0].compare("C") == 0) {	// Capillary
					CAPILLARY_SEGMENT cs;
//					cs.tag = s[1].toInt();
					cs.pos1[0] = s[1].toDouble();
					cs.pos1[1] = s[2].toDouble();
					cs.pos1[2] = s[3].toDouble();
					cs.pos2[0] = s[4].toDouble();
					cs.pos2[1] = s[5].toDouble();
					cs.pos2[2] = s[6].toDouble();
					cs.radius  = s[7].toDouble();
					capillary_list.append(cs);
				} else if (s[0].compare("P") == 0) {
					PIT_POS pp;
					pp.pos[0] = s[1].toInt();
					pp.pos[1] = s[2].toInt();
					pp.pos[2] = s[3].toInt();
					pp.ypit = s[4].toDouble();
					pitpos_list.append(pp);
				} else if (s[0].compare("E") == 0) {
					break;
				} 
			}
		} while (!line.isNull());
	}
	posdata.close();
//	LOG_MSG("VTK: did read posdata");
	if (savepos) {
		delete out;
		vtkdata->close();
		delete vtkdata;
	}

	if (QFile::exists(infileName)) {
		QFile::rename(infileName,"TO_REMOVE");
//		QFile::remove("TO_REMOVE");
//		QFile::remove(infileName);
		LOG_MSG("Did rename");
	}
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
int MyVTK::inTileList(SURFACE_TILE *tile)
{
	for (int i=0; i<tile_list.length(); i++) {
		SURFACE_TILE t = tile_list[i];
		if (t.axis != tile->axis) continue;
		if ((t.pos[0] != tile->pos[0]) || (t.pos[1] != tile->pos[1]) || (t.pos[2] != tile->pos[2])) continue;
		return i;
	}
	return -1;
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::makeBox()
{
	vtkActor *actor;
	double v[3], cpos[3];
	int i, i1, i2, j;
	double boxColor[] = {1,1,1};
	double cnr[8][3];
	int side[8][2];

	cnr[0][0] =      0.5; cnr[0][1] = NBY +0.5; cnr[0][2] =      0.5;
	cnr[1][0] = NX + 0.5; cnr[1][1] = NBY +0.5; cnr[1][2] =      0.5;
	cnr[2][0] = NX + 0.5; cnr[2][1] = NBY +0.5; cnr[2][2] = NZ + 0.5;
	cnr[3][0] =      0.5; cnr[3][1] = NBY +0.5; cnr[3][2] = NZ + 0.5;
	cnr[4][0] =      0.5; cnr[4][1] = NY + 0.5; cnr[4][2] =      0.5;
	cnr[5][0] = NX + 0.5; cnr[5][1] = NY + 0.5; cnr[5][2] =      0.5;
	cnr[6][0] = NX + 0.5; cnr[6][1] = NY + 0.5; cnr[6][2] = NZ + 0.5;
	cnr[7][0] =      0.5; cnr[7][1] = NY + 0.5; cnr[7][2] = NZ + 0.5;
	side[0][0] = 0; side[0][1] = 4;
	side[1][0] = 1; side[1][1] = 5;
	side[2][0] = 2; side[2][1] = 6;
	side[3][0] = 3; side[3][1] = 7;
	side[4][0] = 4; side[4][1] = 5;
	side[5][0] = 5; side[5][1] = 6;
	side[6][0] = 6; side[6][1] = 7;
	side[7][0] = 7; side[7][1] = 4;

	for (i=0; i<8; i++) {
		cnr[i][0] -= 0;
		cnr[i][2] -= 0;
	}

	for(i=0; i<8; i++) {
		i1 = side[i][0];
		i2 = side[i][1];
		actor = vtkActor::New();
		actor->SetMapper(cylMapper);
		actor->GetProperty()->SetColor(boxColor);
		for (j=0; j<3; j++) {
			cpos[j] = (cnr[i1][j] + cnr[i2][j])/2;
			v[j] = cnr[i1][j] - cnr[i2][j];
		}
		double v_mod = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
		double s[] = {1.0, v_mod, 1.0};
		actor->SetScale(s);
		for (j=0; j<3; j++)
			v[j] = v[j]/v_mod;

		double sina = sqrt(v[0]*v[0] + v[2]*v[2]);
		double cosa = v[1];
		double theta = asin(sina)*(180.0/Pi);
		if (cosa < 0)
			theta = 180 - theta;

		actor->SetPosition(cpos);
		actor->RotateWXYZ(theta,v[2],0,-v[0]);
		actor->GetProperty()->SetAmbient(0.5);
		actor->GetProperty()->SetDiffuse(0.2);
		actor->GetProperty()->SetSpecular(0.5);
		ren->AddActor(actor);
		boxActor[i] = actor;
	}
}

//---------------------------------------------------------------------------------------------
// From the list of pits, create a list of bone surface tiles.
// It would really make a lot more sense to start with the full pit depth at each (x,z)
// In this case pos[1] is the last site removed (pit bottom is at pos[1]-1).
//---------------------------------------------------------------------------------------------
void MyVTK::process_tiles()
{
	double pos[3];
	vtkActor *actor;
//	VERT_TILE tile;
//	double v[3];
//	char msg[256];
	unsigned char col[3];
	/*
	double boneColor[] = {0.9,0.9,0.5};
	double pitColor[] = {1.0,0.3,0.1};
	double whiteColor[] = {1,1,1};
	double redColor[] = {1,0,0};
	double greenColor[] = {0,1,0};
	double blueColor[] = {0,0,1};
	double xColor[] = {0.0,1.0,0.0};
	double zColor[] = {0.0,0.0,1.0};
	unsigned char r[3] = {255,0,0};
	unsigned char g[3] = {0,255,0};
	unsigned char b[3] = {0,0,255};
	*/
	unsigned char bone[3] = {250,250,128};

	vtkSmartPointer<vtkCellArray> dummycells;

	// Set up bone surface tile actors.  These are unchanging, except when a new pit is
	// initiated, in which case the position (y) is modified.  This is separate from
	// Tile_Actor_list, which holds all tiles that are not parallel to X-Z.
	// Note that integer site coords (x,y,z) are passed as in (1,NX) etc.
	// The origin of the coord system is actually 0.5 beyond the edge of the lattice boxes.
	// A cell at lattice position (1,1,1) is at the centre of a cube which extends
	// x:0.5-1.5, y:0.5-1.5, z:0.5-1.5.

	if (bone_array == 0) {
		LOG_MSG("allocate new bone_array");
		bone_array = new BONE*[NX+1];
		int npos = NX*NZ;
		points = vtkSmartPointer<vtkPoints>::New();
//		dummycells = vtkSmartPointer<vtkCellArray>::New();
		points->SetNumberOfPoints(npos);         // number of squares
//		dummycells->InsertNextCell(npos);           // size of dummy "verts" array
		colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
//		vtkSmartPointer<vtkDoubleArray> colors = vtkSmartPointer<vtkDoubleArray>::New();
		colors->SetName("colors");
		colors->SetNumberOfComponents(3);
		LOG_MSG("Start loop");
		int x, z, ipos=0;
		for (x = 0; x < NX+1; x++) {
			bone_array[x] = new BONE[NZ+1];
			for (z = 0; z < NZ+1; z++) {
				if (x == 0 || z == 0) {
					bone_array[x][z].actor = 0;
					continue;
				}
				bone_array[x][z].actor = 0;
				pos[0] = x;
				pos[1] = NBY + 0.5;
				pos[2] = z;
				bone_array[x][z].y = pos[1];
				points->SetPoint(ipos, pos[0], pos[1], pos[2]); // square glyph positions set here
//				dummycells->InsertCellPoint(ipos);
				colors->InsertNextTupleValue(bone);
				ipos++;
			}
		}
		LOG_MSG("Made positions, dummycells and colors");
		ginput = vtkSmartPointer<vtkPolyData>::New();
		ginput->SetPoints(points);
//		ginput->SetVerts(dummycells);
		ginput->GetPointData()->SetScalars(colors);
//		vtkSmartPointer<vtkDataArray> c = ginput->GetPointData()->GetScalars();
//		c->SetTuple3( 1000 , 1.0, 0.0, 0.0 );
//		colors->SetTupleValue(1000,zColor);

		if (glypher == 0) {
//			glypher->Delete();
			glypher = vtkSmartPointer<vtkGlyph3D>::New();
		}
		glypher->SetInput(ginput);
		glypher->SetSource(polygonPolyData);
		glypher->ScalingOff();
		glypher->SetColorModeToColorByScalar();
		glypher->Update();

		if (squareMapper == 0) {
//			squareMapper->Delete();
			squareMapper = vtkPolyDataMapper::New();
		}
		squareMapper->SetInputConnection(glypher->GetOutputPort());
		LOG_MSG("did squareMapper->SetInputConnection");
		if (glactor == 0) {
			LOG_MSG("new glactor");
			glactor = vtkActor::New();
		}
		glactor->SetMapper(squareMapper);
		LOG_MSG("did glactor->SetMapper");
	//	glactor->GetProperty()->SetColor(boneColor);
		glactor->GetProperty()->SetAmbient(0.1);	// this washes out the tile colour, -> uniform pink
//		glactor->GetProperty()->SetDiffuse(0.5);	// no effect
//		glactor->GetProperty()->SetSpecular(0.2);
		ren->AddActor(glactor);

		makeBox();
	}

	// tile_list is now used to hold all the vertical faces of the pits.
	// The pitpos_list contains the (x,y,z) of the centre of the bottom of each pit.
	// This is used to modify the y field of each corresponding entry in bone_array[][].
	// When bone_array[][] has been updated to the new geometry, the vertical faces
	// are generated thus:
	// Each (x,z) in pitpos_list is processed in turn.  For each of the 4 neighbour squares,
	// if the neighbour y exceeds the current square's y, a face is created to span the gap,
	// and added to tile_list.  The face actor is made from the basic tile by scaling (in the y
	// direction) and rotation.

	// First need to remove all vertical tile actors and clear the verttile_list
	// (If we are not planning to reuse list entries ...)

	for (int i=0; i<verttile_list.length(); i++) {
		actor = verttile_list[i].actor;
		if (actor != 0) {
			ren->RemoveActor(actor);
			actor->Delete();
		}
	}
	verttile_list.clear();

	//glypher->GetInput();
//	vtkSmartPointer<vtkDataArray> c = ginput->GetPointData()->GetScalars();
	for (int i=0; i<pitpos_list.length(); i++) {
		int x = pitpos_list[i].pos[0];
		int z = pitpos_list[i].pos[2];
		double ypit = pitpos_list[i].ypit;
		bone_array[x][z].y = ypit;
		pos[0] = x;
		pos[1] = ypit;
		pos[2] = z;
// Here is where the glyph position needs to be set
//		bone_array[x][z].actor->SetPosition(pos);
//		bone_array[x][z].actor->GetProperty()->SetColor(pitColor);
		int idx = (x-1)*NZ + z-1;
//		points->SetPoint(idx, pos[0], pos[1], pos[2]); // set point idx to (x,y,z)
//		points->GetData();
		if (ypit == 99) {
			col[1] = 255; col[2] = 255; col[0] = 255;
			colors->SetTupleValue(idx, col);
		} else {
			int redval = 255 - (NBY + 0.5 - ypit)*50;
			if (redval < 0) redval = 0;
			col[1] = 0; col[2] = 0; col[0] = (unsigned char)redval;
			colors->SetTupleValue(idx, col);
		}
		glypher->Modified(); // tell glypher that the data has changed
	}
	return;

	/*
	// The code below, for displaying the sides of the pits, gets slow
	// when there are many pits.
	int nnew = 0;
	for (int i=0; i<pitpos_list.length(); i++) {
		double height;
		int x = pitpos_list[i].pos[0];
//		int ypit = pitpos_list[i].ypit;
		int z = pitpos_list[i].pos[2];
		if (x == 1 || x == NX) continue;
		if (z == 1 || z == NZ) continue;
		height = bone_array[x-1][z].y - bone_array[x][z].y;
		if (height > 0) {
			tile.axis = 'X';
			tile.height = height;
			tile.pos[0] = x - 0.5;
			tile.pos[2] = z;
			tile.pos[1] =(bone_array[x-1][z].y + bone_array[x][z].y)/2;
//			sprintf(msg,"X height: %f pos: %f %f %f",height,pos[0],pos[1],pos[2]);
//			LOG_MSG(msg);
			nnew++;
			tile.actor = vtkActor::New();
			tile.actor->SetMapper(tileMapper);
			verttile_list.append(tile);
//			LOG_MSG("appended tile");
		}
		height = bone_array[x+1][z].y - bone_array[x][z].y;
		if (height > 0) {
			tile.axis = 'X';
			tile.height = height;
			tile.pos[0] = x + 0.5;
			tile.pos[2] = z;
			tile.pos[1] =(bone_array[x+1][z].y + bone_array[x][z].y)/2;
			nnew++;
			tile.actor = vtkActor::New();
			tile.actor->SetMapper(tileMapper);
			verttile_list.append(tile);
		}
		height = bone_array[x][z-1].y - bone_array[x][z].y;
		if (height > 0) {
			tile.axis = 'Z';
			tile.height = height;
			tile.pos[0] = x;
			tile.pos[2] = z - 0.5;
			tile.pos[1] =(bone_array[x][z-1].y + bone_array[x][z].y)/2;
			nnew++;
			tile.actor = vtkActor::New();
			tile.actor->SetMapper(tileMapper);
			verttile_list.append(tile);
		}
		height = bone_array[x][z+1].y - bone_array[x][z].y;
		if (height > 0) {
			tile.axis = 'Z';
			tile.height = height;
			tile.pos[0] = x;
			tile.pos[2] = z + 0.5;
			tile.pos[1] =(bone_array[x][z+1].y + bone_array[x][z].y)/2;
			nnew++;
			tile.actor = vtkActor::New();
			tile.actor->SetMapper(tileMapper);
			verttile_list.append(tile);
		}
	}

	for (int i=0; i<verttile_list.length(); i++) {
		double s[] = {1, 1, 1};
		tile = verttile_list[i];
		if (tile.axis == 'X') {
			s[0] = tile.height;
			v[0] = 0; v[1] = 0; v[2] = 1;
			tile.actor->GetProperty()->SetColor(pitColor);
		} else {
			s[2] = tile.height;
			v[0] = 1; v[1] = 0; v[2] = 0;
			tile.actor->GetProperty()->SetColor(pitColor);
		}
		tile.actor->SetScale(s);
		tile.actor->SetPosition(tile.pos);
		tile.actor->RotateWXYZ(90.0,v[0],v[1],v[2]);
		tile.actor->GetProperty()->SetAmbient(0.5);
		tile.actor->GetProperty()->SetDiffuse(0.2);
		tile.actor->GetProperty()->SetSpecular(0.5);
		ren->AddActor(tile.actor);
	}
	sprintf(msg,"New tiles: %d",nnew);
	LOG_MSG(msg);
	*/
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::init()
{
    T_Actor_list.clear();
	M_Actor_list.clear();
    B_Actor_list.clear();
	C_Actor_list.clear();
	Tile_Actor_list.clear();
	verttile_list.clear();
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::cleanup()
{
	int i;
	vtkActor *actor;
	LOG_MSG("VTK cleanup");
	for (i = 0; i<M_Actor_list.length(); i++) {
		actor = M_Actor_list[i];
		if (actor != 0) {
			ren->RemoveActor(actor);
			actor->Delete();
		}
	}
	for (i = 0; i<C_Actor_list.length(); i++) {
		actor = C_Actor_list[i];
		if (actor != 0) {
			ren->RemoveActor(actor);
			actor->Delete();
		}
	}
	for (i = 0; i<Tile_Actor_list.length(); i++) {
		actor = Tile_Actor_list[i];
		if (actor != 0) {
			ren->RemoveActor(actor);
			actor->Delete();
		}
	}
	if (bone_array != 0) {
		for (int x=0; x<NX+1; x++) {
			for (int z=0; z<NZ+1; z++) {
				actor = bone_array[x][z].actor;
				if (actor != 0) {
					ren->RemoveActor(actor);
					actor->Delete();
				}
			}
			delete bone_array[x];
		}
		delete bone_array;
		bone_array = 0;
		for (i=0; i<8; i++) {
			actor = boxActor[i];
			ren->RemoveActor(actor);
			actor->Delete();
			boxActor[i] = 0;
		}
		LOG_MSG("Deleting glactor");
		ren->RemoveActor(glactor);
		glactor->Delete();
		glactor = 0;
	}
	for (i=0; i<verttile_list.length(); i++) {
		actor = verttile_list[i].actor;
		ren->RemoveActor(actor);
		actor->Delete();
	}

	verttile_list.clear();
	T_Actor_list.clear();
	M_Actor_list.clear();
    B_Actor_list.clear();
	C_Actor_list.clear();
	Tile_Actor_list.clear();
	first_VTK = true;
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::renderCells()
{
//	LOG_MSG("VTK: renderCells");

	process_Mcells();
	process_capillaries();
	process_tiles();
	process_osteoclasts();

	if (first_VTK) {
		LOG_MSG("Initializing the renderer");
		ren->ResetCamera();
//		iren->Render();
	}
	iren->Render();
	first_VTK = false;
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::process_Tcells()
{
	int i, tag;
	double r, g, b, genfac;
//	double activated = 9.0;
	double TC_MAX_GEN = 15;
	CELL_POS cp;
	vtkActor *actor;
	double TCColor[] = {0.0, 0.0, 1.0};

    int na = T_Actor_list.length();
    int np = TCpos_list.length();
    int n = na;
	for (i=0; i<np; i++) {
        cp = TCpos_list[i];
        tag = cp.tag;
        n = max(tag+1,n);
	}
    bool *active;
	active = new bool[n];
	for (i=0; i<n; i++)
		active[i] = false;
	for (i=0; i<np; i++) {
        cp = TCpos_list[i];
        tag = cp.tag;
        active[tag] = true;
		if (tag >= na) {   // need to add actor, and possibly fill gaps
			if (tag > na) {
                for (int j=na; j<tag; j++)	//j in range(na,tag):
                    T_Actor_list.append(0);
			}
			actor = vtkActor::New();
            actor->SetMapper(TcellMapper);
            actor->GetProperty()->SetColor(TCColor);
//            actor->GetProperty()->SetColor(0.0,0.0,1.0);	// TCColor
            ren->AddActor(actor);
            T_Actor_list.append(actor);
            na = tag + 1;
//			sprintf(msg,"added T_actor: %d  %p  %p",tag,ren,actor);
//			LOG_MSG(msg);
		}
		if (cp.state == -1) {	// non-cognate
			r = 0.5; g = 0.5; b = 0.5;
		} else if (cp.state == 0) {
			r = 0; g = 0; b = 1;
		} else if (cp.state <= TC_MAX_GEN) {
			genfac = (cp.state-1)/(TC_MAX_GEN-1);		// 0 - 1
			b = genfac*0.4;
			g = 1 - b;
			r = 0;
		} else {
			r = 1.0; g = 0.6; b = 0.0;
		}

        actor = T_Actor_list[tag];
        actor->GetProperty()->SetColor(r, g, b);
        actor->SetPosition(cp.x, cp.y, cp.z);
		if (actor != 0) 
			actor->SetPosition(cp.x, cp.y, cp.z);
		else {
			sprintf(msg,"T_actor = 0: %d",tag);
			LOG_MSG(msg);
			exit(1);
		}
	}

	for (int k=0; k<na; k++) {	// k in range(0,na):
		if (T_Actor_list[k] != 0 && !active[k]) {     // need to remove actor from list
            actor = T_Actor_list[k];
            ren->RemoveActor(actor);
			actor->Delete();
			T_Actor_list[k] = 0;
		}
	}
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::process_Mcells()
{
	int i, tag;
	CELL_POS cp;
	vtkActor *actor;
	double r=0, g=0, b=0;
	int nf = 0;
	int nnew = 0;

	int na = M_Actor_list.length();
	int np = MCpos_list.length();
//	sprintf(msg,"process_Mcells: na, np: %d %d",na,np);
//	LOG_MSG(msg);
    int n = na;
	for (i=0; i<np; i++) {
		cp = MCpos_list[i];
        tag = cp.tag;
        n = max(tag+1,n);
	}
    bool *active;
	active = new bool[n];
	for (i=0; i<n; i++)
		active[i] = false;
	int nactive = 0;
	for (i=0; i<np; i++) {
		cp = MCpos_list[i];
        tag = cp.tag;
        active[tag] = true;
		nactive++;
		bool newMC = false;
		if (tag >= na) {   // need to add actor, and possibly fill gaps
			if (tag > na) {
                for (int j=na; j<tag; j++)	//j in range(na,tag):
					M_Actor_list.append(0);
			}
			nnew++;
			actor = vtkActor::New();
			actor->SetMapper(McellMapper);
//			actor->GetProperty()->SetColor(MCColor);
//            actor->GetProperty()->SetColor(0.7,0.2,0.3);	// DCColor

            ren->AddActor(actor);
			M_Actor_list.append(actor);
            na = tag + 1;
			newMC = true;
//			sprintf(msg,"added M_actor: %d",tag);
//			LOG_MSG(msg);
		} else {
			actor = M_Actor_list[tag];
			if (actor == 0) {
				nnew++;
				actor = vtkActor::New();
				actor->SetMapper(McellMapper);
				ren->AddActor(actor);
				M_Actor_list[tag] = actor;
				newMC = true;
			}
		}
		if (actor != 0) {
			actor->SetPosition(cp.x, cp.y, cp.z);
			if (cp.state == 0) {			// MOTILE
				r = 0; g = 0.5; b = 0.5;
			} else if (cp.state == 1) {		// CROSSING
				r = 0; g = 1; b = 0;
			} else if (cp.state < 100) {	// FUSING
				r = (154 + cp.state)/255.; g = r; b = 0;	//b = (100 - cp.state)/255.;
			} else if (cp.state == 100) {						// FUSED
				nf++;
				r = 1; g = 1; b = 0;
			}
			if (tag%100 == 0) {				// stain some cells to observe motility
				r = 1; g = 1; b = 1;
			}
			actor->GetProperty()->SetColor(r, g, b);
		} else {
//			sprintf(msg,"M_actor = 0: %d  %d",tag,cp.state);
//			LOG_MSG(msg);
//			exit(1);
		}
	}

	for (int k=0; k<na; k++) {	// k in range(0,na):
		if (M_Actor_list[k] != 0 && !active[k]) {     // need to remove actor from list
			actor = M_Actor_list[k];
            ren->RemoveActor(actor);
			actor->Delete();
			M_Actor_list[k] = 0;
		}
	}
//	sprintf(msg,"nnew, nactive: %d %d",nnew,nactive);
//	LOG_MSG(msg);
}


//---------------------------------------------------------------------------------------------
// A cylinder is created orientated along the y-axis, i.e. along b = (0,1,0)
// To change the orientation to the vector v, we first create a vector r
// normal to both b and v: r = bxv, this will be the axis of rotation.
// We now need to rotate the cylinder by theta about r, where theta is the angle between
// b and v, i.e. sin(theta) = |r|/(|b||v|) = |r|/|v|
// We can now use actor.RotateWXYZ(theta,r[0],r[1],r[2]) where theta is in degrees
// What is bxv when b = (0,1,0) and v = (v0,v1,v2)?
// r = [v[2],0,-v[0]]
//---------------------------------------------------------------------------------------------
void MyVTK::process_capillaries()
{
	int i, j;
	CAPILLARY_SEGMENT cs;
	vtkActor *actor;
	double cpos[3], v[3];
	double capillaryColor[] = {0.5,0.0,0.2};

	int na = C_Actor_list.length();
	int np = capillary_list.length();

	// First remove all old capillaries (strictly speaking we should remove only those not in the new list)

	for (int k=0; k<na; k++) {
		actor = C_Actor_list[k];
		ren->RemoveActor(actor);
		actor->Delete();
		C_Actor_list[k] = 0;
	}

	C_Actor_list.clear();

	for (i=0; i<np; i++) {
		cs = capillary_list[i];
		actor = vtkActor::New();
		actor->SetMapper(capillaryMapper);
		actor->GetProperty()->SetColor(capillaryColor);
		for (j=0; j<3; j++) {
			cpos[j] = (cs.pos1[j] + cs.pos2[j])/2;
			v[j] = cs.pos1[j] - cs.pos2[j];
		}
		double v_mod = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
		double s[] = {cs.radius, v_mod, cs.radius};
		actor->SetScale(s);
		for (j=0; j<3; j++)
			v[j] = v[j]/v_mod;

		double sina = sqrt(v[0]*v[0] + v[2]*v[2]);
		double cosa = v[1];
		double theta = asin(sina)*(180.0/Pi);
		if (cosa < 0)
			theta = 180 - theta;

		actor->SetPosition(cpos);
		actor->RotateWXYZ(theta,v[2],0,-v[0]);
		ren->AddActor(actor);
		C_Actor_list.append(actor);
	}
}

//---------------------------------------------------------------------------------------------
void MyVTK::process_osteoclasts()
{
	int iclast, k;
//	CAPILLARY_SEGMENT cs;
	CLAST_POS oc;
	vtkActor *actor;
	double theta;
	double clastColor[] = {0.9,0.9,0.0};

	int na = OC_Actor_list.length();
	int np = clastpos_list.length();

	// First remove all old clasts (strictly speaking we should remove only those not in the new list)

	for (k=0; k<na; k++) {
		actor = OC_Actor_list[k];
		ren->RemoveActor(actor);
		actor->Delete();
		OC_Actor_list[k] = 0;
	}

	OC_Actor_list.clear();

	for (iclast=0; iclast<np; iclast++) {
		oc = clastpos_list[iclast];
		actor = vtkActor::New();
		actor->SetMapper(clastMapper);
		actor->GetProperty()->SetColor(clastColor);
		theta = (180./Pi)*atan2(oc.dir[2],oc.dir[0]);
//		sprintf(msg,"x,z,theta: %d %f %f %f",iclast,oc.dir[0],oc.dir[2],-theta);
//		LOG_MSG(msg);
		actor->RotateWXYZ(-theta,0,1,0);

		actor->SetScale(oc.size);
		actor->SetPosition(oc.pos);
		ren->AddActor(actor);
		OC_Actor_list.append(actor);
	}
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool MyVTK::startPlayer(QString posfile, QTimer *theTimer, bool save)
{
//	casename = acasename;
	save_image = save;
//	posfile = casename + ".pos";
	LOG_QMSG(posfile);
	timer = theTimer;
	playerData = new QFile(posfile);
	if (!playerData->open(QFile::ReadOnly)) {
		LOG_MSG("Open failure on VTK file");
		return false;
	}
	playerStream = new QTextStream(playerData);
	if (!first_VTK) {
		cleanup();
	}
	playing = true;
	paused = false;

	if (save_image) {
		w2i = vtkWindowToImageFilter::New();
		w2i->SetInput(renWin);	//the render window
//		writer = vtkSmartPointer<vtkPNGWriter>::New();
		writer = vtkSmartPointer<vtkJPEGWriter>::New();

//		castFilter = vtkSmartPointer<vtkImageCast>::New();
//		castFilter->SetOutputScalarTypeToUnsignedChar ();
//		castFilter->SetInputConnection(w2i->GetOutputPort());
//		castFilter->Update();

		writer->SetInputConnection(w2i->GetOutputPort()); 
		framenum = 0;
		LOG_MSG("set up writer");
	}

	LOG_MSG("playing");
	return true;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool MyVTK::nextFrame()
{
	LOG_MSG("VTK: nextFrame");
	if (!playing)
		return false;
	if (paused)
		return true;
	if (playerStream->atEnd()) {
		LOG_MSG("nextFrame: no more data");
		stop();
		return false;
	}
	TCpos_list.clear();
	MCpos_list.clear();
//	bondpos_list.clear();
	capillary_list.clear();
	tile_list.clear();
	int k = 0;
	QString line;
	do {
		line = playerStream->readLine();
		if (line.length() > 0) {
			k++;
			QStringList s = line.split(" ",QString::SkipEmptyParts);
			if (s[0].compare("T") == 0) {
				CELL_POS cp;
				cp.tag = s[1].toInt();
				cp.x = s[2].toInt();
				cp.y = s[3].toInt();
				cp.z = s[4].toInt();
//				cp.diameter = s[5].toDouble();
//				cp.state = s[6].toDouble();
				TCpos_list.append(cp);
			} else if (s[0].compare("M") == 0) {
				CELL_POS cp;
				cp.tag = s[1].toInt();
				cp.x = s[2].toInt();
				cp.y = s[3].toInt();
				cp.z = s[4].toInt();
//				cp.diameter = s[5].toDouble();
				cp.state = s[5].toInt();
				MCpos_list.append(cp);
			} else if (s[0].compare("B") == 0) {
				BOND_POS cp;
				cp.TCtag = s[1].toInt();
				cp.DCtag = s[2].toInt();
				bondpos_list.append(cp);
			} else if (s[0].compare("E") == 0) {
				break;
			}
		}
	} while (true);

	bool redo = false;
	if (first_VTK) {
		redo = true;
	}
//    renderCells(redo,false);
	renderCells();
	char numstr[5];
	sprintf(numstr,"%04d",framenum);
	if (save_image) {
		w2i->Modified();	//importante 
		writer->SetFileName((casename + numstr + ".jpg").toStdString().c_str()); 
		writer->Write(); 
	}
	framenum++;
	return true;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::saveSnapshot(QString fileName, QString imgType)
{
//	vtkWindowToImageFilter *w2img = vtkWindowToImageFilter::New();
	w2img->SetInput(renWin);
	if (imgType.compare("png") == 0) {
		vtkSmartPointer<vtkPNGWriter> pngwriter = vtkPNGWriter::New();
		pngwriter->SetInputConnection(w2img->GetOutputPort()); 
		w2img->Modified();
		pngwriter->SetFileName((fileName.toStdString()).c_str()); 
		pngwriter->Write();
//		pngwriter->Delete();	// Note: using vtkSmartPointer, delete is not necessary.
	} else if (imgType.compare("jpg") == 0) {
		vtkJPEGWriter *jpgwriter = vtkJPEGWriter::New();
		jpgwriter->SetInputConnection(w2img->GetOutputPort()); 
		w2img->Modified();
		jpgwriter->SetFileName((fileName.toStdString()).c_str()); 
		jpgwriter->Write();
//		jpgwriter->Delete();
	} else if (imgType.compare("tif") == 0) {
		vtkTIFFWriter *tifwriter = vtkTIFFWriter::New();
		tifwriter->SetInputConnection(w2img->GetOutputPort()); 
		w2img->Modified();
		tifwriter->SetFileName((fileName.toStdString()).c_str()); 
		tifwriter->Write();
//		tifwriter->Delete();
	} else if (imgType.compare("bmp") == 0) {
		vtkBMPWriter *bmpwriter = vtkBMPWriter::New();
		bmpwriter->SetInputConnection(w2img->GetOutputPort()); 
		w2img->Modified();
		bmpwriter->SetFileName((fileName.toStdString()).c_str()); 
		bmpwriter->Write();
//		bmpwriter->Delete();
	}
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::pause()
{
	paused = true;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::playon()
{
	paused = false;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::stop()
{
	if (save_image) {
		writer->Delete();
		w2i->Delete();
	}
	delete playerStream;
	playerData->close();
	delete playerData;
	timer->stop();
	playing = false;
	paused = false;
}

