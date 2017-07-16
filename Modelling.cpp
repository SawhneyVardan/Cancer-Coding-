//
//  Modelling.cpp
//  
//
//  Created by Vardan Sawhney on 2017-01-06.
//
//

#include "vtkCubeSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkPolyData.h"
#include "vtkGlyph3DMapper.h"

vtkSmartPointer<vtkRenderer> renderer;
vtkSmartPointer<vtkRenderWindow> renderWindow;
vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
vtkSmartPointer<vtkPoints> points;

void initializeRendering() {
    renderer = vtkSmartPointer<vtkRenderer>::New(); //creating new renderer
    renderer->SetBackground(1., 1., 1.); //setting background color to white
    
    renderWindow = vtkSmartPointer<vtkRenderWindow>::New(); //creating renderer window
    renderWindow->SetSize(1000, 1000); //setting size of the window
    renderWindow->AddRenderer(renderer); //adding previously created renderer
    
    renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New(); //creating render
    renderWindowInteractor->SetRenderWindow(renderWindow); //adding previously created renderer window
    
    points = vtkSmartPointer<vtkPoints>::New();         //vector with points
}

void addCells(vector<Cell>::iterator begin, vector<Cell>::iterator end) {
    vector<Cell>::iterator cell = begin;
    for (cell; cell<end; cell++)
        points->InsertNextPoint(cell->x,cell->y,cell->z);
}

void visualizeLesion(vector<Cell>::iterator begin, vector<Cell>::iterator end) {
    
    initializeRendering();
    
    addCells(begin, end);
    
    // Combine into a polydata
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    
    vtkSmartPointer<vtkCubeSource> cubeSource = vtkSmartPointer<vtkCubeSource>::New();//you can use sphere if you want
    
    vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
    glyph3D->SetSourceConnection(cubeSource->GetOutputPort());
    glyph3D->ScalingOff();
    glyph3D->SetInputData(polydata);
    glyph3D->Update();
    
    // Create a mapper and actor
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(glyph3D->GetOutputPort());
    
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    
    // Add the actor to the scene
    renderer->AddActor(actor);
    
    // Render
    renderWindow->Render();
    
    vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New(); //like paraview
    
    renderWindowInteractor->SetInteractorStyle( style );
    
    renderWindowInteractor->Start();
}

void addCells(vector<Cell>::iterator begin, vector<Cell>::iterator end) {
    vector<Cell>::iterator cell = begin;
    for (cell; cell<end; cell++) {         if (cell->z < 250) {             points->InsertNextPoint(cell->x,cell->y,cell->z);
    } else {
        points->InsertNextPoint(cell->x,cell->y,cell->z+50);
    }
    }
}

#include "vtkCubeSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkPolyData.h"
#include "vtkGlyph3DMapper.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkLookupTable.h"

vtkSmartPointer<vtkRenderer> renderer;
vtkSmartPointer<vtkRenderWindow> renderWindow;
vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;

vtkSmartPointer<vtkPoints> points;

vtkSmartPointer<vtkFloatArray> scalar;


void initializeRendering() {
    renderer = vtkSmartPointer<vtkRenderer>::New(); //creating new renderer
    renderer->SetBackground(1., 1., 1.); //setting background color to white
    
    renderWindow = vtkSmartPointer<vtkRenderWindow>::New(); //creating renderer window
    renderWindow->SetSize(1000, 1000); //setting size of the window
    renderWindow->AddRenderer(renderer); //adding previously created renderer
    
    renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New(); //creating render
    renderWindowInteractor->SetRenderWindow(renderWindow); //adding previously created renderer window
    
    points = vtkSmartPointer<vtkPoints>::New();         //vector with points
    scalar = vtkSmartPointer<vtkFloatArray>::New(); //vector with scalar values of the points
    scalar->SetNumberOfComponents(1); //one value per point
}

void addCells(vector<Cell>::iterator begin, vector<Cell>::iterator end) {
    vector<Cell>::iterator cell = begin;
    for (cell; cell<end; cell++) {         if (cell->z < 250) {             points->InsertNextPoint(cell->x,cell->y,cell->z);
    } else {
        points->InsertNextPoint(cell->x,cell->y,cell->z+50);
    }
        
        scalar->InsertNextValue((double)cell->nSpots/26.);
    }
}

void visualizeLesion(vector<Cell>::iterator begin, vector<Cell>::iterator end) {
    
    initializeRendering();
    
    addCells(begin, end);
    
    //creating new color table
    vtkSmartPointer<vtkLookupTable> colorLookupTable = vtkSmartPointer<vtkLookupTable>::New();
    colorLookupTable->SetNumberOfColors(256);
    colorLookupTable->SetHueRange( 0.667, 0.0);
    colorLookupTable->Build();
    
    // Combine into a polydata
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    polydata->GetPointData()->SetScalars(scalar);
    
    vtkSmartPointer<vtkCubeSource> cubeSource = vtkSmartPointer<vtkCubeSource>::New();//you can use sphere if you want
    
    vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
    glyph3D->SetColorModeToColorByScalar();
    glyph3D->SetSourceConnection(cubeSource->GetOutputPort());
    glyph3D->ScalingOff();
    glyph3D->SetInputData(polydata);
    glyph3D->Update();
    
    // Create a mapper and actor
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(glyph3D->GetOutputPort());
    mapper->SetLookupTable(colorLookupTable);
    
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    
    // Add the actor to the scene
    renderer->AddActor(actor);
    
    // Render
    renderWindow->Render();
    
    vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New(); //like paraview
    
    renderWindowInteractor->SetInteractorStyle( style );
    
    renderWindowInteractor->Start();
}
