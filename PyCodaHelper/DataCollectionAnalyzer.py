import numpy as np
import os
import vtk

class DataCollectionAnalyzer:
    def __init__(self, baseFolder):
        self.__baseFolder = baseFolder
        self.__vrtxToPrimVrtxMap = dict()
        self.__elemToPrimElemMap = dict()
        
    def Help(self):
        raise Exception("Coming Soon")
    
    def CheckSymmetry(self,
                      ptDataArrNames,
                      cellDataArrNames,
                      symPlanes,
                      absTh = 1e-4,
                      serOrMpi = 'serial',
                      DCFolder="ParaView",
                      cycles=['./']):#, 'Cycle000001']):
        # checks
        if(serOrMpi not in ['serial']):#, 'mpi', 'both']
            raise Exception("serOrMpi has to be among serial/mpi/both")
        
        for cycle in cycles:
            for root, dirs, files in os.walk(self.__baseFolder+'/'+DCFolder+'/'+cycle):
                for name in files:
                    fileName = os.path.join(root, name)
                    if(serOrMpi == 'serial') and (fileName.endswith('.vtu')):
                        print(fileName)
                        self.__CheckDataSymmetry(fileName, 
                                                 ptDataArrNames, 
                                                 cellDataArrNames, symPlanes, absTh)
    def __CreateSymMaps(self, fileName, symPlanes):

        # Unstructured Grid reader
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(fileName)
        reader.Update()

        # Remove duplicate points
        clnData = vtk.vtkStaticCleanUnstructuredGrid()
        clnData.SetInputData(reader.GetOutput())
        clnData.Update()

        # Unstructured Grid
        ug = clnData.GetOutput()

        # Some reusable parameters
        numPts = ug.GetNumberOfPoints()
        numCells = ug.GetNumberOfCells()
        
        # create point set from primary vertices
        primVrtxPts = vtk.vtkPoints()
        ctr = 0
        zero = 1e-12
        vrtxMapVtkIdToMfemId = dict()
        for i in range(numPts):
            vertex = np.array(ug.GetPoint(i))
            # check if the given point is in primary domain
            isPrimaryVertex = True
            for symPlane in symPlanes:
                func_val = symPlane.EvalPlaneEquation(vertex)
                if(func_val>zero): #the given vertex is not a primary vertex
                    isPrimaryVertex = False
                    break

            if(isPrimaryVertex):
                primVrtxPts.InsertNextPoint(vertex)
                vrtxMapVtkIdToMfemId[ctr] = i
                ctr += 1

        # create point set from primary element centroids
        primElemCentPts = vtk.vtkPoints()
        ctr = 0
        elemMapVtkIdToMfemId = dict()
        for i in range(numCells):
            cell = ug.GetCell(i)

            # Compute the center of the cell
            cellPts = cell.GetPoints()
            numCellPts = cellPts.GetNumberOfPoints()
            centroid = np.zeros(3)
            for j in range(numCellPts):
                centroid += np.array(cellPts.GetPoint(j))
            centroid /= numCellPts

            # check if the given centroid is in primary domain
            isPrimaryElement = True
            for symPlane in symPlanes:
                func_val = symPlane.EvalPlaneEquation(centroid)
                if(func_val>zero): #the given centroid is not in primary domain
                    isPrimaryElement = False
                    break

            if(isPrimaryElement):
                primElemCentPts.InsertNextPoint(centroid)
                elemMapVtkIdToMfemId[ctr] = i
                ctr += 1

        # Create KD Tree from primary domain vertices
        kdTreePrimVrts = vtk.vtkKdTree()
        kdTreePrimVrts.BuildLocatorFromPoints(primVrtxPts)

        # Create KD Tree from primary element centroids
        kdTreePrimElems = vtk.vtkKdTree()
        kdTreePrimElems.BuildLocatorFromPoints(primElemCentPts)

        # Create map from vertex to primary vertex id
        dist = vtk.reference(100000000.0)
        result = vtk.vtkIdList()
        vrtxToPrimVrtxMap = dict()
        for i in range(numPts):
            pt = ug.GetPoint(i)
            
            # Check if pt is in the primary half of the full mesh
            # If it's not in the primary half then reflect it (sometimes it requires
            # multiple reflections) to get it into the primary half
            for symPlane in symPlanes:
                func_val = symPlane.EvalPlaneEquation(np.array(pt));
                if(func_val > zero): # if val is greater than ZERO then it's not in the primary half
                    pt = symPlane.ReflectPoint(np.array(pt));          
            
            p = kdTreePrimVrts.FindClosestPoint(pt, dist)
            #print("My point : {} | reference point : {} | {}".format(ug.GetPoint(i), primVrtxPts.GetPoint(p), dist )) 
            self.__vrtxToPrimVrtxMap[i] = vrtxMapVtkIdToMfemId[p]

        # Create map from element to primary element i
        elemToPrimElemMap = dict()
        for i in range(numCells):
            cell = ug.GetCell(i)

            # Compute the center of the cell
            cellPts = cell.GetPoints()
            numCellPts = cellPts.GetNumberOfPoints()
            centroid = np.zeros(3)
            for j in range(numCellPts):
                centroid += np.array(cellPts.GetPoint(j))
            centroid /= numCellPts
            
            # Check if pt is in the primary half of the full mesh
            # If it's not in the primary half then reflect it (sometimes it requires
            # multiple reflections) to get it into the primary half
            for symPlane in symPlanes:
                func_val = symPlane.EvalPlaneEquation(centroid);
                if(func_val > zero): # if val is greater than ZERO then it's not in the primary half
                    centroid = symPlane.ReflectPoint(np.array(centroid));    

            # Find the primary element
            p = kdTreePrimElems.FindClosestPoint(centroid, dist)
            self.__elemToPrimElemMap[i] = elemMapVtkIdToMfemId[p]
                        
    def __CheckDataSymmetry(self, fileName, ptDataArrNames, cellDataArrNames, symPlanes, absTh):
        # Create Symmetry maps
        self.__CreateSymMaps(fileName, symPlanes)
        
        # Unstructured Grid reader
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(fileName)
        reader.Update()

        # Remove duplicate points
        clnData = vtk.vtkStaticCleanUnstructuredGrid()
        clnData.SetInputData(reader.GetOutput())
        clnData.Update()

        # Unstructured Grid
        ug = clnData.GetOutput()

        # Some reusable parameters
        numPts = ug.GetNumberOfPoints()
        numCells = ug.GetNumberOfCells()
        
        # Check fror symmetry in pointData
        ptData = ug.GetPointData()
        numViolateVrtx = 0
        for arrName in ptDataArrNames:
            arr = ptData.GetArray(arrName)
            for i in range(numPts):
                myVal = arr.GetTuple1(i)
                refVal = arr.GetTuple1(self.__vrtxToPrimVrtxMap[i])
                absErr = abs(myVal-refVal)
                if(absErr >= absTh):
                    print("my coord : {} | reference coord : {} | {}-array my values : {} | reference values :{}  Difference : {} \n".format(
                        ug.GetPoint(i), ug.GetPoint(self.__vrtxToPrimVrtxMap[i]), arrName, myVal, refVal, absErr))
                    numViolateVrtx += 1
        
        # Check fror symmetry in cellData
        cellData = ug.GetCellData()
        numViolateElem = 0
        for arrName in cellDataArrNames:
            arr = cellData.GetArray(arrName)
            for i in range(numCells):
                myVal = arr.GetTuple1(i)
                refVal = arr.GetTuple1(self.__elemToPrimElemMap[i])
                absErr = abs(myVal-refVal)
                if (absErr >= absTh): 
                    bounds = ug.GetCell(i).GetBounds()
                    print("cell bounds: {}, {}, {}, {}, {}, {} | {}-array my values : {} | reference values :{}  Difference : {} \n".format(
                        bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5], arrName, myVal, refVal, absErr))
                    numViolateElem += 1
                    
        if (numViolateVrtx > 0) or (numViolateElem > 0):
            return False
        return True