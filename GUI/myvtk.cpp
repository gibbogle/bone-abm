// myvtk.cpp

#ifdef _WIN32
#include "windows.h"
#endif
#include "myvtk.h"
#include "log.h"
#include "transfer.h"

LOG_USE();

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
MyVTK::MyVTK(QWidget *page)
{
	zoomlevel = 0.7;
	double backgroundColor[] = {0.0,0.0,0.0};
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
	vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
	points->InsertNextPoint(-0.5, 0.0, -0.5);
	points->InsertNextPoint( 0.5, 0.0, -0.5);
	points->InsertNextPoint( 0.5, 0.0,  0.5);
	points->InsertNextPoint(-0.5, 0.0,  0.5);
	// a unit square in the XY plane, at y=0

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

	//Create a mapper and actor
	tileMapper = vtkPolyDataMapper::New();
	tileMapper->SetInput(polygonPolyData);
	tileMapper->ScalarVisibilityOff();

	vtkSmartPointer<vtkTextSource> textSource = vtkSmartPointer<vtkTextSource>::New();
	textSource->SetText("Hello");
	textSource->SetForegroundColor(1.0, 1.0, 1.0);
	textSource->BackingOff();
	textSource->Update();

	// Create a mapper and actor
	textMapper = vtkPolyDataMapper::New();
	textMapper->SetInputConnection(textSource->GetOutputPort());

	//legendScaleActor = vtkSmartPointer<vtkLegendScaleActor>::New();

	// Create image filter for save Snapshot()
	w2img = vtkWindowToImageFilter::New();
//	pngwriter = vtkSmartPointer<vtkPNGWriter>::New();
//	jpgwriter = vtkSmartPointer<vtkJPEGWriter>::New();

	bone_array = 0;
	first_VTK = true;
	DCmotion = false;
	DCfade = true;
	playing = false;
	paused = false;
	Pi = 4*atan(1.0);

	ren->GetActiveCamera()->Zoom(zoomlevel);		// try zooming OUT
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
MyVTK::~MyVTK()
{
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::get_cell_positions()
{
	LOG_MSG("VTK: get_cell_positions");
	MCpos_list.clear();
	capillary_list.clear();
	pitpos_list.clear();
	for (int i=0; i<nmono_list; i++) {
		int j = 5*i;
		CELL_POS cp;
		cp.tag = mono_list[j];
		cp.x = mono_list[j+1];
		cp.y = mono_list[j+2];
		cp.z = mono_list[j+3];
		cp.state = mono_list[j+4];
		MCpos_list.append(cp);
	}
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
	double pos[3];
	vtkActor *actor;
	VERT_TILE tile;
	double s[3], v[3], cpos[3];
	int i, i1, i2, j;
	char msg[256];
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
	VERT_TILE tile;
	double s[3], v[3];
	char msg[256];
	double boneColor[] = {0.9,0.9,0.5};
	double pitColor[] = {0.9,0.4,0.2};
	double xColor[] = {0.0,1.0,0.0};
	double zColor[] = {0.0,0.0,1.0};

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
		for (int x = 0; x < NX+1; x++) {
			bone_array[x] = new BONE[NZ+1];
			for (int z = 0; z < NZ+1; z++) {
				if (x == 0 || z == 0) {
					bone_array[x][z].actor = 0;
					continue;
				}
				actor = vtkActor::New();
				actor->SetMapper(tileMapper);
				actor->GetProperty()->SetColor(boneColor);
				actor->GetProperty()->SetAmbient(0.5);
				actor->GetProperty()->SetDiffuse(0.2);
				actor->GetProperty()->SetSpecular(0.5);
				pos[0] = x;
				pos[1] = NBY + 0.5;
				pos[2] = z;
				actor->SetPosition(pos);
				ren->AddActor(actor);
				bone_array[x][z].y = pos[1];
				bone_array[x][z].actor = actor;
			}
		}
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

	for (int i=0; i<pitpos_list.length(); i++) {
		int x = pitpos_list[i].pos[0];
		int y = pitpos_list[i].pos[1];
		int z = pitpos_list[i].pos[2];
		double ypit = pitpos_list[i].ypit;
		bone_array[x][z].y = ypit;
		pos[0] = x;
		pos[1] = ypit;
		pos[2] = z;
		bone_array[x][z].actor->SetPosition(pos);
		bone_array[x][z].actor->GetProperty()->SetColor(pitColor);
	}
	for (int i=0; i<pitpos_list.length(); i++) {
		double height;
		int x = pitpos_list[i].pos[0];
		int ypit = pitpos_list[i].ypit;
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
			tile.actor = vtkActor::New();
			tile.actor->SetMapper(tileMapper);
			verttile_list.append(tile);
		}
		height = bone_array[x+1][z].y - bone_array[x][z].y;
		if (height > 0) {
			tile.axis = 'X';
			tile.height = height;
			tile.pos[0] = x + 0.5;
			tile.pos[2] = z;
			tile.pos[1] =(bone_array[x+1][z].y + bone_array[x][z].y)/2;
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
}

//---------------------------------------------------------------------------------------------
// From the list of pits, create a list of bone surface tiles.
// It would really make a lot more sense to start with the full pit depth at each (x,z)
// In this case pos[1] is the last site removed (pit bottom is at pos[1]-1).
//---------------------------------------------------------------------------------------------
void MyVTK::oldcreateTileList()
{
	int j, y, ypit;
	SURFACE_TILE tile;
	char msg[256];

	tile_list.clear();
	y = NBY;

	for (int x=1; x<=NX; x++) {
		for (int z=1; z<=NZ; z++) {
			tile.pos[0] = x;
			tile.pos[1] = y;
			tile.pos[2] = z;
			tile.axis = 'Y';
			tile.pit = false;
			tile_list.append(tile);
		}
	}

	for (int i=0; i<pitpos_list.length(); i++) {
		for (j=0; j<3; j++)
			tile.pos[j] = pitpos_list[i].pos[j];
		ypit = tile.pos[1];

		// First find the Y-tile corresponding to y = NBY, and put pos[1] = ypit-1
		tile.axis = 'Y';
		tile.pos[1] = NBY;
		j = inTileList(&tile);
		if (j >= 0) {
			tile_list[j].pos[1] = ypit-1;
			tile_list[j].pit = true;
		} else {
			sprintf(msg,"Error: createTileList: tile not found: pos: %d %d %d axis: %c\n",tile.pos[0],tile.pos[1],tile.pos[2],tile.axis);
			LOG_MSG(msg);
			exit(1);
		}
		// Now either add or remove tiles at each vertical face
		for (y=ypit; y<=NBY; y++) {
			tile.pit = true;
			tile.pos[1] = y;
			tile.axis = 'X';
			j = inTileList(&tile);
			if (j >= 0)
				tile_list.removeAt(j);
			else
				tile_list.append(tile);
			tile.pos[0] = tile.pos[0]-1;
			j = inTileList(&tile);
			if (j >= 0)
				tile_list.removeAt(j);
			else
				tile_list.append(tile);
			tile.pos[0] = tile.pos[0]+1;

			tile.axis = 'Z';
			j = inTileList(&tile);
			if (j >= 0)
				tile_list.removeAt(j);
			else
				tile_list.append(tile);
			tile.pos[2] = tile.pos[2]-1;
			j = inTileList(&tile);
			if (j >= 0)
				tile_list.removeAt(j);
			else
				tile_list.append(tile);
			tile.pos[2] = tile.pos[2]+1;
		}
	}
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

	if (first_VTK) {
		LOG_MSG("Initializing the renderer");
		ren->ResetCamera();
	}
	for(int i=0; i<1; i++) {
		iren->Render();
	}
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

	int na = M_Actor_list.length();
	int np = MCpos_list.length();
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

	for (i=0; i<np; i++) {
		cp = MCpos_list[i];
        tag = cp.tag;
        active[tag] = true;
		bool newMC = false;
		if (tag >= na) {   // need to add actor, and possibly fill gaps
			if (tag > na) {
                for (int j=na; j<tag; j++)	//j in range(na,tag):
					M_Actor_list.append(0);
			}
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
		}
		if (actor != 0) {
			actor->SetPosition(cp.x, cp.y, cp.z);
			if (cp.state == 0) {			// MOTILE
				r = 0; g = 0.5; b = 0.5;
			} else if (cp.state == 1) {		// CROSSING
				r = 0; g = 1; b = 0;
			} else if (cp.state < 101) {	// FUSING
				r = (154 + cp.state)/255.; g = 1; b = (101 - cp.state)/255.;
			} else {						// FUSED
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
}

/*
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
void MyVTK::process_bonds()
{
	int i, j;
	BOND_POS bp;
	vtkActor *actor, *T_actor, *D_actor;
	double bpos[3], v[3];
	double Pi = 3.15159;
	double *tcpos, *dcpos;
	double bondColor[] = {0.5,0.0,0.0};

    int na = B_Actor_list.length();
    int np = bondpos_list.length();

    // First remove all old bonds (strictly speaking we should remove only those not in the new list)

	for (int k=0; k<na; k++) {
		actor = B_Actor_list[k];
		ren->RemoveActor(actor);
		actor->Delete();
	}

    B_Actor_list.clear();    

	for (i=0; i<np; i++) {
        bp = bondpos_list[i];
		actor = vtkActor::New();
        actor->SetMapper(bondMapper);
		actor->GetProperty()->SetColor(bondColor);
//        actor->GetProperty()->SetColor(0.7,0.2,0.3);	// bondColor
		T_actor = T_Actor_list[bp.TCtag];
		if (T_actor != 0)
	        tcpos = T_actor->GetPosition();
		else {
			sprintf(msg,"T_actor = 0 in bond: %d %d",i,bp.TCtag);
			LOG_MSG(msg);
			exit(1);
		}
		D_actor = D_Actor_list[bp.DCtag];
		if (D_actor != 0)
	        dcpos = D_actor->GetPosition();
		else {
			sprintf(msg,"D_actor = 0 in bond: %d %d",i,bp.DCtag);
			LOG_MSG(msg);
			exit(1);
		}
	
		for (j=0; j<3; j++) {
            bpos[j] = (tcpos[j] + dcpos[j])/2;
            v[j] = tcpos[j] - dcpos[j];
		}
        double v_mod = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
		double s[] = {1, v_mod, 1};
        actor->SetScale(s);
        for (j=0; j<3; j++)
            v[j] = v[j]/v_mod;
            
        double sina = sqrt(v[0]*v[0] + v[2]*v[2]);
        double cosa = v[1];
        double theta = asin(sina)*(180.0/Pi);
		if (cosa < 0) 
            theta = 180 - theta;
		
        actor->SetPosition(bpos);
        actor->RotateWXYZ(theta,v[2],0,-v[0]);
        ren->AddActor(actor);
        B_Actor_list.append(actor);
	}
}
*/

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
	double capillaryColor[] = {0.5,0.0,0.0};

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
void MyVTK::oldprocess_tiles()
{
	int i;
	SURFACE_TILE tile;
	vtkActor *actor;
	double pos[3], v[3];
	double boneColor[] = {0.9,0.9,0.5};
	double pitColor[] = {0.9,0.4,0.2};
	char msg[256];

	oldcreateTileList();
	int na = Tile_Actor_list.length();
	int np = tile_list.length();

	sprintf(msg,"na, np: %d %d",na,np);
	LOG_MSG(msg);

	// First remove all old tiles
	for (int k=0; k<na; k++) {
		actor = Tile_Actor_list[k];
		ren->RemoveActor(actor);
		actor->Delete();
		Tile_Actor_list[k] = 0;
	}
	Tile_Actor_list.clear();

	for (i=0; i<np; i++) {
		tile = tile_list[i];
		actor = vtkActor::New();
		actor->SetMapper(tileMapper);
		if (tile.pit)
			actor->GetProperty()->SetColor(pitColor);
		else
			actor->GetProperty()->SetColor(boneColor);

		// experimenting ...
		actor->GetProperty()->SetAmbient(0.5);
		actor->GetProperty()->SetDiffuse(0.2);
		actor->GetProperty()->SetSpecular(0.5);

		pos[0] = tile.pos[0];
		pos[1] = tile.pos[1];
		pos[2] = tile.pos[2];
		if (tile.axis == 'X') {
			pos[0] = tile.pos[0] + 0.5;
			v[0] = 0; v[1] = 0; v[2] = 1;
		} else if (tile.axis == 'Y') {
			pos[1] = tile.pos[1] + 0.5;
		} else if (tile.axis == 'Z') {
			pos[2] = tile.pos[2] + 0.5;
			v[0] = 1; v[1] = 0; v[2] = 0;
		}
		actor->SetPosition(pos);
		if (tile.axis == 'X' || tile.axis == 'Z') {
			actor->RotateWXYZ(90.0,v[0],v[1],v[2]);
		}
		ren->AddActor(actor);
		Tile_Actor_list.append(actor);
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

