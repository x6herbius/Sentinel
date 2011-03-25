//
//  MapWindowController.m
//  TrenchBroom
//
//  Created by Kristian Duske on 15.03.10.
//  Copyright 2010 __MyCompanyName__. All rights reserved.
//

#import "MapWindowController.h"
#import "MapView3D.h"
#import "TextureView.h"
#import "MapDocument.h"
#import "Entity.h"
#import "Brush.h"
#import "Face.h"
#import "Camera.h"
#import "MapDocument.h"
#import "WadLoader.h"
#import "Wad.h"
#import "TextureManager.h"
#import "InputManager.h"
#import "VBOBuffer.h"
#import "Octree.h"
#import "Picker.h"
#import "SelectionManager.h"
#import "GLFontManager.h"
#import "InspectorController.h"
#import "Options.h"
#import "Grid.h"
#import "Ray3D.h"
#import "Vector3f.h"
#import "Vector3i.h"
#import "PrefabManager.h"
#import "PrefabNameSheetController.h"
#import "Prefab.h"
#import "MapWriter.h"
#import "CameraAnimation.h"
#import "TrackingManager.h"

static NSString* CameraDefaults = @"Camera";
static NSString* CameraDefaultsFov = @"Field Of Vision";
static NSString* CameraDefaultsNear = @"Near Clipping Plane";
static NSString* CameraDefaultsFar = @"Far Clipping Plane";

@implementation MapWindowController

- (NSUndoManager *)windowWillReturnUndoManager:(NSWindow *)window {
    return [[self document] undoManager];
}

- (void)windowDidBecomeKey:(NSNotification *)notification {
    InspectorController* inspector = [InspectorController sharedInspector];
    [inspector setMapWindowController:self];
}

- (void)userDefaultsChanged:(NSNotification *)notification {
    NSDictionary* cameraDefaults = [[NSUserDefaults standardUserDefaults] dictionaryForKey:CameraDefaults];
    if (cameraDefaults == nil)
        return;
    
    float fov = [[cameraDefaults objectForKey:CameraDefaultsFov] floatValue];
    float near = [[cameraDefaults objectForKey:CameraDefaultsNear] floatValue];
    float far = [[cameraDefaults objectForKey:CameraDefaultsFar] floatValue];
    
    [camera setFieldOfVision:fov];
    [camera setNearClippingPlane:near];
    [camera setFarClippingPlane:far];
}

- (void)windowDidLoad {
    GLResources* glResources = [[self document] glResources];
    NSOpenGLContext* sharedContext = [glResources openGLContext];
    NSOpenGLContext* sharingContext = [[NSOpenGLContext alloc] initWithFormat:[view3D pixelFormat] shareContext:sharedContext];
    [view3D setOpenGLContext:sharingContext];
    [sharingContext release];
    
    options = [[Options alloc] init];
    camera = [[Camera alloc] init];
    [self userDefaultsChanged:nil];
    
    selectionManager = [[SelectionManager alloc] init];
    trackingManager = [[TrackingManager alloc] initWithWindowController:self];
    inputManager = [[InputManager alloc] initWithWindowController:self];
    
    [view3D setup];
    
    InspectorController* inspector = [InspectorController sharedInspector];
    [inspector setMapWindowController:self];

    NSNotificationCenter* center = [NSNotificationCenter defaultCenter];
    [center addObserver:self selector:@selector(windowDidBecomeKey:) name:NSWindowDidBecomeKeyNotification object:[self window]];
    [center addObserver:self selector:@selector(userDefaultsChanged:) name:NSUserDefaultsDidChangeNotification object:[NSUserDefaults standardUserDefaults]];
    
    [[self window] setAcceptsMouseMovedEvents:YES];
    [[self window] makeKeyAndOrderFront:nil];
}

- (Camera *)camera {
    return camera;
}

- (SelectionManager *)selectionManager {
    return selectionManager;
}

- (InputManager *)inputManager {
    return inputManager;
}

- (TrackingManager *)trackingManager {
    return trackingManager;
}

- (Options *)options {
    return options;
}

- (BOOL)validateMenuItem:(NSMenuItem *)menuItem {
    SEL action = [menuItem action];
    if (action == @selector(clearSelection:)) {
        return [selectionManager hasSelection];
    } else if (action == @selector(copySelection:)) {
        return [selectionManager hasSelectedEntities] || [selectionManager hasSelectedBrushes];
    } else if (action == @selector(cutSelection:)) {
        return [selectionManager hasSelectedEntities] || [selectionManager hasSelectedBrushes];
    } else if (action == @selector(pasteClipboard:)) {
        return NO;
    } else if (action == @selector(deleteSelection:)) {
        return [selectionManager hasSelectedEntities] || [selectionManager hasSelectedBrushes];
    } else if (action == @selector(moveTextureLeft:)) {
        return [selectionManager hasSelectedBrushes] || [selectionManager hasSelectedFaces];
    } else if (action == @selector(moveTextureLeft:)) {
        return [selectionManager hasSelectedBrushes] || [selectionManager hasSelectedFaces];
    } else if (action == @selector(moveTextureRight:)) {
        return [selectionManager hasSelectedBrushes] || [selectionManager hasSelectedFaces];
    } else if (action == @selector(moveTextureUp:)) {
        return [selectionManager hasSelectedBrushes] || [selectionManager hasSelectedFaces];
    } else if (action == @selector(moveTextureDown:)) {
        return [selectionManager hasSelectedBrushes] || [selectionManager hasSelectedFaces];
    } else if (action == @selector(stretchTextureHorizontally:)) {
        return [selectionManager hasSelectedBrushes] || [selectionManager hasSelectedFaces];
    } else if (action == @selector(shrinkTextureHorizontally:)) {
        return [selectionManager hasSelectedBrushes] || [selectionManager hasSelectedFaces];
    } else if (action == @selector(stretchTextureVertically:)) {
        return [selectionManager hasSelectedBrushes] || [selectionManager hasSelectedFaces];
    } else if (action == @selector(shrinkTextureVertically:)) {
        return [selectionManager hasSelectedBrushes] || [selectionManager hasSelectedFaces];
    } else if (action == @selector(rotateTextureLeft:)) {
        return [selectionManager hasSelectedBrushes] || [selectionManager hasSelectedFaces];
    } else if (action == @selector(rotateTextureRight:)) {
        return [selectionManager hasSelectedBrushes] || [selectionManager hasSelectedFaces];
    } else if (action == @selector(duplicateSelection:)) {
        return [selectionManager hasSelectedBrushes];
    } else if (action == @selector(createPrefabFromSelection:)) {
        return [selectionManager hasSelectedBrushes];
    } else if (action == @selector(showInspector:)) {
        return YES;
    } else if (action == @selector(switchToXYView:)) {
        return YES;
    } else if (action == @selector(switchToInvertedXYView:)) {
        return YES;
    } else if (action == @selector(switchToXZView:)) {
        return YES;
    } else if (action == @selector(switchToInvertedXZView:)) {
        return YES;
    } else if (action == @selector(switchToYZView:)) {
        return YES;
    } else if (action == @selector(switchToInvertedYZView:)) {
        return YES;
    } else if (action == @selector(isolateSelection:)) {
        return YES;
    } else if (action == @selector(toggleProjection:)) {
        return YES;
    }

    return NO;
}

- (IBAction)showInspector:(id)sender {
    InspectorController* inspector = [InspectorController sharedInspector];
    [[inspector window] makeKeyAndOrderFront:nil];
}

- (IBAction)toggleGrid:(id)sender {
    [[options grid] toggleDraw];
}

- (IBAction)toggleSnap:(id)sender {
    [[options grid] toggleSnap];
}

- (IBAction)setGridSize:(id)sender {
    [[options grid] setSize:[sender tag]];
}

- (IBAction)isolateSelection:(id)sender {
    EIsolationMode isolationMode = [options isolationMode];
    [options setIsolationMode:(isolationMode + 1) % 3];
}

- (IBAction)toggleProjection:(id)sender {
    ECameraMode cameraMode = [camera mode];
    [camera setMode:(cameraMode + 1) % 2];
}

- (IBAction)switchToXYView:(id)sender {
    Vector3f* center = [selectionManager selectionCenter];
    if (center == nil) {
        center = [[Vector3f alloc] initWithFloatVector:[camera direction]];
        [center scale:256];
        [center add:[camera position]];
    }
    
    Vector3f* diff = [[Vector3f alloc] initWithFloatVector:center];
    [diff sub:[camera position]];
    
    Vector3f* position = [[Vector3f alloc] initWithFloatVector:center];
    [position setZ:[position z] + [diff length]];
    
    CameraAnimation* animation = [[CameraAnimation alloc] initWithCamera:camera targetPosition:position targetDirection:[Vector3f zAxisNeg] targetUp:[Vector3f yAxisPos] duration:0.5];
    [animation startAnimation];
    
    [diff release];
    [position release];
}

- (IBAction)switchToXZView:(id)sender {
    Vector3f* center = [selectionManager selectionCenter];
    if (center == nil) {
        center = [[Vector3f alloc] initWithFloatVector:[camera direction]];
        [center scale:256];
        [center add:[camera position]];
    }
    
    Vector3f* diff = [[Vector3f alloc] initWithFloatVector:center];
    [diff sub:[camera position]];
    
    Vector3f* position = [[Vector3f alloc] initWithFloatVector:center];
    [position setY:[position y] - [diff length]];
    
    CameraAnimation* animation = [[CameraAnimation alloc] initWithCamera:camera targetPosition:position targetDirection:[Vector3f yAxisPos] targetUp:[Vector3f zAxisPos] duration:0.5];
    [animation startAnimation];

    [diff release];
    [position release];
}

- (IBAction)switchToYZView:(id)sender {
    Vector3f* center = [selectionManager selectionCenter];
    if (center == nil) {
        center = [[Vector3f alloc] initWithFloatVector:[camera direction]];
        [center scale:256];
        [center add:[camera position]];
    }
    
    Vector3f* diff = [[Vector3f alloc] initWithFloatVector:center];
    [diff sub:[camera position]];
    
    Vector3f* position = [[Vector3f alloc] initWithFloatVector:center];
    [position setX:[position x] + [diff length]];
    
    CameraAnimation* animation = [[CameraAnimation alloc] initWithCamera:camera targetPosition:position targetDirection:[Vector3f xAxisNeg] targetUp:[Vector3f zAxisPos] duration:0.5];
    [animation startAnimation];
    
    [diff release];
    [position release];
}

- (IBAction)clearSelection:(id)sender {
    [selectionManager removeAll];
}

- (IBAction)copySelection:(id)sender {}
- (IBAction)cutSelection:(id)sender {}
- (IBAction)pasteClipboard:(id)sender {}

- (IBAction)deleteSelection:(id)sender {
    NSUndoManager* undoManager = [[self document] undoManager];
    [undoManager beginUndoGrouping];

    NSSet* deletedBrushes = [[NSSet alloc] initWithSet:[selectionManager selectedBrushes]];

    NSEnumerator* brushEn = [deletedBrushes objectEnumerator];
    id <Brush> brush;
    while ((brush = [brushEn nextObject]))
        [[self document] deleteBrush:brush];
    
    [selectionManager removeAll];
    [deletedBrushes release];
    
    [undoManager endUndoGrouping];
    [undoManager setActionName:@"Delete Selection"];
}

- (IBAction)moveTextureLeft:(id)sender {
    NSUndoManager* undoManager = [[self document] undoManager];
    [undoManager beginUndoGrouping];
    
    int d = ![[options grid] snap] ^ ([NSEvent modifierFlags] & NSAlternateKeyMask) != 0 ? 1 : [[options grid] size];

    NSEnumerator* faceEn = [[selectionManager selectedFaces] objectEnumerator];
    id <Face> face;
    while ((face = [faceEn nextObject]))
        [[self document] translateFaceOffset:face xDelta:-d yDelta:0];

    [undoManager endUndoGrouping];
    [undoManager setActionName:@"Move Texture"];
}

- (IBAction)moveTextureRight:(id)sender {
    NSUndoManager* undoManager = [[self document] undoManager];
    [undoManager beginUndoGrouping];

    int d = ![[options grid] snap] ^ ([NSEvent modifierFlags] & NSAlternateKeyMask) != 0 ? 1 : [[options grid] size];
    
    NSEnumerator* faceEn = [[selectionManager selectedFaces] objectEnumerator];
    id <Face> face;
    while ((face = [faceEn nextObject]))
        [[self document] translateFaceOffset:face xDelta:d yDelta:0];
    
    [undoManager endUndoGrouping];
    [undoManager setActionName:@"Move Texture"];
}

- (IBAction)moveTextureUp:(id)sender {
    NSUndoManager* undoManager = [[self document] undoManager];
    [undoManager beginUndoGrouping];

    int d = ![[options grid] snap] ^ ([NSEvent modifierFlags] & NSAlternateKeyMask) != 0 ? 1 : [[options grid] size];
    
    NSEnumerator* faceEn = [[selectionManager selectedFaces] objectEnumerator];
    id <Face> face;
    while ((face = [faceEn nextObject]))
        [[self document] translateFaceOffset:face xDelta:0 yDelta:d];
    
    [undoManager endUndoGrouping];
    [undoManager setActionName:@"Move Texture"];
}

- (IBAction)moveTextureDown:(id)sender {
    NSUndoManager* undoManager = [[self document] undoManager];
    [undoManager beginUndoGrouping];

    int d = ![[options grid] snap] ^ ([NSEvent modifierFlags] & NSAlternateKeyMask) != 0 ? 1 : [[options grid] size];
    
    NSEnumerator* faceEn = [[selectionManager selectedFaces] objectEnumerator];
    id <Face> face;
    while ((face = [faceEn nextObject]))
        [[self document] translateFaceOffset:face xDelta:0 yDelta:-d];
    
    [undoManager endUndoGrouping];
    [undoManager setActionName:@"Move Texture"];
}

- (IBAction)stretchTextureHorizontally:(id)sender {
    NSUndoManager* undoManager = [[self document] undoManager];
    [undoManager beginUndoGrouping];

    NSEnumerator* faceEn = [[selectionManager selectedFaces] objectEnumerator];
    id <Face> face;
    while ((face = [faceEn nextObject]))
        [[self document] setFace:face xScale:[face xScale] + 0.1f];
    
    [undoManager endUndoGrouping];
    [undoManager setActionName:@"Stretch Texture Horizontally"];
}

- (IBAction)shrinkTextureHorizontally:(id)sender {
    NSUndoManager* undoManager = [[self document] undoManager];
    [undoManager beginUndoGrouping];
    
    NSEnumerator* faceEn = [[selectionManager selectedFaces] objectEnumerator];
    id <Face> face;
    while ((face = [faceEn nextObject]))
        [[self document] setFace:face xScale:[face xScale] - 0.1f];
    
    [undoManager endUndoGrouping];
    [undoManager setActionName:@"Shrink Texture Horizontally"];
}

- (IBAction)stretchTextureVertically:(id)sender {
    NSUndoManager* undoManager = [[self document] undoManager];
    [undoManager beginUndoGrouping];
    
    NSEnumerator* faceEn = [[selectionManager selectedFaces] objectEnumerator];
    id <Face> face;
    while ((face = [faceEn nextObject]))
        [[self document] setFace:face yScale:[face yScale] + 0.1f];
    
    [undoManager endUndoGrouping];
    [undoManager setActionName:@"Stretch Texture Vertically"];
}

- (IBAction)shrinkTextureVertically:(id)sender {
    NSUndoManager* undoManager = [[self document] undoManager];
    [undoManager beginUndoGrouping];
    
    NSEnumerator* faceEn = [[selectionManager selectedFaces] objectEnumerator];
    id <Face> face;
    while ((face = [faceEn nextObject]))
        [[self document] setFace:face yScale:[face yScale] - 0.1f];
    
    [undoManager endUndoGrouping];
    [undoManager setActionName:@"Shrink Texture Vertically"];
}

- (IBAction)rotateTextureLeft:(id)sender {
    NSUndoManager* undoManager = [[self document] undoManager];
    [undoManager beginUndoGrouping];
    
    int d = ![[options grid] snap] ^ ([NSEvent modifierFlags] & NSAlternateKeyMask) != 0 ? 1 : 15;
    
    NSEnumerator* faceEn = [[selectionManager selectedFaces] objectEnumerator];
    id <Face> face;
    while ((face = [faceEn nextObject]))
        [[self document] setFace:face rotation:[face rotation] - d];
    
    [undoManager endUndoGrouping];
    [undoManager setActionName:@"Rotate Texture Left"];
}

- (IBAction)rotateTextureRight:(id)sender {
    NSUndoManager* undoManager = [[self document] undoManager];
    [undoManager beginUndoGrouping];
    
    int d = ![[options grid] snap] ^ ([NSEvent modifierFlags] & NSAlternateKeyMask) != 0 ? 1 : 15;
    
    NSEnumerator* faceEn = [[selectionManager selectedFaces] objectEnumerator];
    id <Face> face;
    while ((face = [faceEn nextObject]))
        [[self document] setFace:face rotation:[face rotation] + d];
    
    [undoManager endUndoGrouping];
    [undoManager setActionName:@"Rotate Texture Right"];
}

- (IBAction)duplicateSelection:(id)sender {
    NSUndoManager* undoManager = [[self document] undoManager];
    [undoManager beginUndoGrouping];
    
    id <Entity> worldspawn = [[self document] worldspawn];
    NSMutableSet* newBrushes = [[NSMutableSet alloc] init];

    NSEnumerator* brushEn = [[selectionManager selectedBrushes] objectEnumerator];
    id <Brush> brush;
    while ((brush = [brushEn nextObject])) {
        id <Brush> newBrush = [[self document] createBrushInEntity:worldspawn fromTemplate:brush];
        [[self document] translateBrush:newBrush xDelta:[[options grid] size] yDelta:[[options grid] size] zDelta:[[options grid] size]];
        [newBrushes addObject:newBrush];
    }
    
    [[undoManager prepareWithInvocationTarget:selectionManager] addBrushes:[NSSet setWithSet:[selectionManager selectedBrushes]]];
    [[undoManager prepareWithInvocationTarget:selectionManager] removeAll];
    
    [selectionManager removeAll];
    [selectionManager addBrushes:newBrushes];
    [newBrushes release];
    
    [undoManager endUndoGrouping];
    [undoManager setActionName:@"Duplicate Selection"];
}

- (void)prefabNameSheetDidEnd:(NSWindow *)sheet returnCode:(NSInteger)returnCode contextInfo:(void *)contextInfo {
    PrefabNameSheetController* pns = [sheet windowController];
    if (returnCode == NSOKButton) {
        NSString* prefabName = [pns prefabName];
        NSString* prefabGroupName = [pns prefabGroup];
        
        PrefabManager* prefabManager = [PrefabManager sharedPrefabManager];
        id <PrefabGroup> prefabGroup = [prefabManager prefabGroupWithName:prefabGroupName create:YES];
        [prefabManager createPrefabFromBrushTemplates:[selectionManager selectedBrushes] name:prefabName group:prefabGroup];
    }
    [pns release];
}

- (IBAction)createPrefabFromSelection:(id)sender {
    PrefabNameSheetController* pns = [[PrefabNameSheetController alloc] init];
    NSWindow* prefabNameSheet = [pns window];
    
    NSApplication* app = [NSApplication sharedApplication];
    [app beginSheet:prefabNameSheet modalForWindow:[self window] modalDelegate:self didEndSelector:@selector(prefabNameSheetDidEnd:returnCode:contextInfo:) contextInfo:nil];
    
}

- (void)insertPrefab:(id <Prefab>)prefab {
    [selectionManager removeAll];
    
    NSUndoManager* undoManager = [[self document] undoManager];
    [undoManager beginUndoGrouping];
    
    Vector3f* insertPos = [[Vector3f alloc] initWithFloatVector:[camera direction]];
    [insertPos scale:256];
    [insertPos add:[camera position]];
    [[options grid] snapToGrid:insertPos];
    
    Vector3f* offset = [[options grid] gridOffsetOf:[prefab center]];
    [insertPos add:offset];

    Vector3f* dist = [[Vector3f alloc] initWithFloatVector:insertPos];
    [dist sub:[prefab center]];
    
    MapDocument* map = [self document];
    
    NSEnumerator* entityEn = [[prefab entities] objectEnumerator];
    id <Entity> prefabEntity;
    while ((prefabEntity = [entityEn nextObject])) {
        id <Entity> mapEntity;
        if ([prefabEntity isWorldspawn]) {
            mapEntity = [map worldspawn];
            if (mapEntity == nil)
                mapEntity = [map createEntityWithProperties:[NSDictionary dictionaryWithObject:@"worldspawn" forKey:@"classname"]];
        } else {
            mapEntity = [map createEntityWithProperties:[prefabEntity properties]];
            [selectionManager addEntity:mapEntity];
        }
        
        NSEnumerator* prefabBrushEn = [[prefabEntity brushes] objectEnumerator];
        id <Brush> prefabBrush;
        while ((prefabBrush = [prefabBrushEn nextObject])) {
            id <Brush> mapBrush = [map createBrushInEntity:mapEntity fromTemplate:prefabBrush];
            [map translateBrush:mapBrush xDelta:roundf([dist x]) yDelta:roundf([dist y]) zDelta:roundf([dist z])];
            [selectionManager addBrush:mapBrush];
        }
    }
    
    [dist release];
    [insertPos release];

    [undoManager endUndoGrouping];
    [undoManager setActionName:[NSString stringWithFormat:@"Insert Prefab '%@'", [prefab name]]];
}

- (void)dealloc {
    [[NSNotificationCenter defaultCenter] removeObserver:self];
    [options release];
    [trackingManager release];
    [selectionManager release];
    [inputManager release];
    [camera release];
    [super dealloc];
}

@end
