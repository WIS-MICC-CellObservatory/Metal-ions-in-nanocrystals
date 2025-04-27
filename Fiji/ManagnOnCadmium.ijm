
//String(label="Process Mode", choices=("singleFile", "wholeFolder", "AllSubFolders"), style="list") gProcessMode
var gProcessMode = "singleFile"; 

#@ File(label="Segmentation image file",value="A:\\ranabu\\CdS-Mn Project\\Dataset 1\\without pre-filter\\REA-CdS-Mn-NP 1159 140 kx SI HAADF-BF_original Mn-net.tif", persist=true, description="for segmenting the particles and analyzing their base distribution") iCellposFile
#@ File(label="Particles image file",value="A:\\ranabu\\CdS-Mn Project\\Dataset 1\\without pre-filter\\REA-CdS-Mn-NP 1159 140 kx SI HAADF-BF_original Cd-net.tif", persist=true, description="for segmenting the particles and analyzing their base distribution") iParticlesFile
#@ File (label="Mangan image file",value="A:\\ranabu\\CdS-Mn Project\\Dataset 1\\without pre-filter\\REA-CdS-Mn-NP 1159 140 kx SI HAADF-BF_original Mn-net.tif", persist=true, description="For analyzing Mn distribution") iMnFile
#@ String(label="Mangan minimal intensity",value="0.5", persist=true, description="below this number Mn values will be ignored'") iMnMinThreshold
#@ File(label="Segmentation image file",value="A:\\ranabu\\CdS-Mn Project\\Dataset 1\\without pre-filter\\REA-CdS-Mn-NP 1159 140 kx SI HAADF-BF_original Mn-net.tif", persist=true, description="for segmenting the particles and analyzing their base distribution") iCellposFile
// String(label="Particle atom minimal intensity",value="0", persist=true, description="below this number Mn values will be ignored'") iparticleMinThreshold
#@ String(label="Results directory",value="Results", persist=true, description="sub directory for results'") iResultsDir
#@ File(label="Cellpose, Enviorenment",value="D:\\ProgramData\\Anaconda3\\envs\\Cellpose", style="directory", persist=true, description="Cellpose env.") iCellposeEnv
#@ String(label="Cellpose, Cell diameter",value="10", persist=true, description="as set in training") iParticleCellposeCellDiameter
#@ boolean(label="Cellpose, Use previous run",value=true, persist=true, description="for quicker runs, use previous labeling of image, if exists") iUseCellposePrevRun

/*
 * Given 3 files: File to segment
 */

//----Macro parameters-----------
var pMacroName = "ManganOnCadmium";
var pMacroVersion = "1.0.1";

// global file variables
var gFileFullPath = "uninitialized1";
var gFileNameNoExt = "uninitialized2";
var gResultsSubFolder = "uninitialized3";
var gImagesResultsSubFolder = "uninitialized4";
var gMainDirectory = "uninitialized5";
var gSubFolderName = "";
var gSaveRunDir = "SaveRun";


//----- global variables-----------
var gCadmiumCellposeModel = "cyto2";
var gCellposeExtRunSubdir = "/Segmentation/";
var gCellposeExtRunFileSuffix = "uninitialized6";

var	gCompositeTable = "CompositeResults";
var	gAllCompositeTable = "allCompositeTable";
var gAllCompositeResults = 0; // the comulative number of rows in allCompositeTable

var width, height, channels, slices, frames;
var unit,pixelWidth, pixelHeight;


//-----debug variables--------
var gDebugFlag = false;
var gBatchModeFlag = false;
//------------Main------------
Initialization();
if(LoopFiles())
	print("Macro ended successfully");
else
	print("Macro failed");
CleanUp(true);
waitForUser("=================== Done ! ===================");

	
function ProcessFile(directory) 
{
	gImagesResultsSubFolder = gResultsSubFolder;
	setBatchMode(gBatchModeFlag);
	initVariables();
	
	open(iMnFile);
	gMnImageId = getImageID(); 
		
	open(iCellposFile);
	gCellposeImageId = getImageID(); 
	rename("segmentationImage");
	
	open(iParticlesFile);
	giParticlesImageId = getImageID(); 
	rename("segmentationImage");
	
	getDimensions(width, height, channels, slices, frames);
	getPixelSize(unit,pixelWidth, pixelHeight);	

	gFileNameNoExt = File.getNameWithoutExtension(iCellposFile);	
	gImagesResultsSubFolder = gResultsSubFolder;	
	File.makeDirectory(gImagesResultsSubFolder);
	File.makeDirectory(gImagesResultsSubFolder+"/"+gSaveRunDir);
	
	// get Particles segmentation using cellpose
	gLabeledParticlesImageId = RunCellposeModel(gCellposeImageId, gCadmiumCellposeModel);
	run("glasbey on dark");
	saveAs("Tiff",gImagesResultsSubFolder+"/Segmentation.tif");

	// generate roi for each particle and calculate its area
	run("Label image to ROIs");
	roiManager("measure");

	// get distance image
	gDistanceImageId = GetDistanceImage(gLabeledParticlesImageId);
	saveAs("Tiff",gImagesResultsSubFolder+"/Distance.tif");

	// get the max distance
	selectImage(gDistanceImageId);
	getRawStatistics(nPixels, mean, min, max, std, histogram);
	maxParticleDepth = Math.ceil(max)+1;	
	numParticlesPerDepth = newArray(maxParticleDepth);
	
	// 1. go over all the particles
	// 2. set in ParticeDepth array, the maximal depth of each particle
	// 3. create a depth table for each new maximal depth (row for the number of particles in each 
	// 4. Update in numParticlesPerDepth array the number of particles with each maximal depth found
	n = roiManager("count");
	ParticleDepth = newArray(n);
	selectImage(gDistanceImageId);
	for(i=0;i<n;i++)
	{
		roiManager("select", i);
		getRawStatistics(nPixels, mean, min, max, std, histogram);
		depth = Math.ceil(max);
		ParticleDepth[i] = depth;	
		if(numParticlesPerDepth[depth] <= 0)
		{
			tableName = gCompositeTable+depth;
			Table.create(tableName);	
			InitRows(tableName,1,depth+1);	
		}
		numParticlesPerDepth[depth]++;
	}
	setBatchMode(true);
	Table.create(gCompositeTable);
	//Table.create(gAllCompositeTable);
	lastLabel = 0;
	// go over all pixels in the labeled image
	// if the label is non zero:
	// 1. get its depth from the distance image
	// 2. increment the number of particles in the relevant depth table
	// 3. in the 
	
	
	for(y=0;y<height;y++)
	{
		for(x=0;x<width;x++)
		{
			selectImage(gLabeledParticlesImageId);
			pId = getPixel(x, y);
			if(pId > 0)
			{
				if(pId>lastLabel)
				{
					// create rows in table for each distance
					InitRows(gCompositeTable,pId,maxParticleDepth);
					//if(pId == 1)
					//	InitRows(gAllCompositeTable,pId,maxParticleDepth);					
					lastLabel = pId;
				}						
				selectImage(gDistanceImageId);
				depth = Math.ceil(getPixel(x, y));
				aggCompositeTable = gCompositeTable+ParticleDepth[pId-1];
				IncTableValue(gCompositeTable,"Pixel count",(pId-1)*maxParticleDepth+depth,1);
				IncTableValue(aggCompositeTable,"Pixel count",depth,1);
				
				selectImage(giParticlesImageId);
				pIntensity = getPixel(x, y);
				IncTableValue(gCompositeTable,"Pixels sum intensities",(pId-1)*maxParticleDepth+depth,pIntensity);
				IncTableValue(aggCompositeTable,"Pixels sum intensities",depth,pIntensity);
				
				selectImage(gMnImageId);
				mnIntensity = getPixel(x, y);
				if(mnIntensity > iMnMinThreshold)
				{
					IncTableValue(gCompositeTable,"Mn count",(pId-1)*maxParticleDepth+depth,1);
					IncTableValue(aggCompositeTable,"Mn count",depth,1);
	
					IncTableValue(gCompositeTable,"Mns sum intensities",(pId-1)*maxParticleDepth+depth,mnIntensity);
					IncTableValue(aggCompositeTable,"Mns sum intensities",depth,mnIntensity);	
				}
			}
		}
	}
	setBatchMode(gBatchModeFlag);	
	Table.save(gImagesResultsSubFolder+"/paritcal.csv",gCompositeTable);
	for(depth=0;depth<maxParticleDepth;depth++)
	{
		aggCompositeTable = gCompositeTable+depth;
		if(numParticlesPerDepth[depth] > 0)
			Table.save(gImagesResultsSubFolder+"/"+aggCompositeTable+"_"+numParticlesPerDepth[depth]+".csv",gCompositeTable+depth);	
	}
	Table.save(gImagesResultsSubFolder+"/ParticleArea.csv","Results");	
	return true;
}
function IncTableValue(table,ColName,tableRow,incValue)
{
	prevValue = Table.get(ColName, tableRow ,table);
	Table.set(ColName, tableRow ,prevValue+incValue,table);		
}

function InitRows(table,pId,numRowsPerParticle)
{
	for(i=0;i<numRowsPerParticle;i++)
	{						
		Table.set("Particle Id", (pId-1)*numRowsPerParticle+i,pId,table);
		Table.set("Distance from Edge", (pId-1)*numRowsPerParticle+i,i,table);
		Table.set("Pixel count", (pId-1)*numRowsPerParticle+i,0,table);				
		Table.set("Pixels sum intensities", (pId-1)*numRowsPerParticle+i,0,table);				
		Table.set("Mn count", (pId-1)*numRowsPerParticle+i,0,table);				
		Table.set("Mns sum intensities", (pId-1)*numRowsPerParticle+i,0,table);				
	}
}
function GetDistanceImage(labelImageId)
{
	ErodeLabels(labelImageId);
	setThreshold(1, 65535, "raw");
	run("Convert to Mask");
	run("Invert");
	run("Distance Transform 3D");
	return getImageID();
}


function initVariables()
{
	//gFirstClusterRoi = true;
}

function GetCellsRois()
{
	gManualRoi = false;
	manualRoiPath = gImagesResultsSubFolder+"/"+gCellsRois+"_Manual.zip";
	if(File.exists(manualRoiPath))
	{
		print("Warning: Using user generated cells rois");
		gManualRoi = true;
		roiManager("open", manualRoiPath);
		return;
	}

	gCellsLabelsImageId = RunCellposeModel(gMitoZimageId, gMitoCellposeModel);
	gNumCells = GenerateROIsFromLabelImage(gCellsLabelsImageId,"Cell",0);
	if(gNumCells <= 0)
	{
		title = getTitle();
		print("WARNING!!!: Cellpose did not identify any cell/object in " + title);
	}
	StoreROIs(gImagesResultsSubFolder,gCellsRois);	
}


function GetLDsRois()
{
	gLDsLabelsImageId = RunCellposeModel(gLDZimageId, gLDCellposeModel);
	gNumLDs = GenerateROIsFromLabelImage(gLDsLabelsImageId,"",0);
	if(gNumLDs <= 0)
	{
		print("WARNING!!!: Cellpose did not identify any LD in " + getTitle());
	}
	StoreROIs(gImagesResultsSubFolder,gLDRois);	
}

function ProcessCell(roiId)
{
	tmpROIFullPath = gImagesResultsSubFolder+"/"+gTempROIFile;
	if(!matches(gProcessMode, "singleFile"))
		Table.set("File Name",gAllCompositeResults+roiId,gFileFullPath,gAllCompositeTable);
		
	SaveROIs(tmpROIFullPath);

	SetCompositeTables("Cell ID",roiId,"Cell_"+(roiId+1));
	SetCompositeTables("Manual Cell roi",roiId,gManualRoi);

	
	//check LD clustering in cell
	AnalyzeLDs(roiId);
	ClearRoiManager();
	roiManager("Open", tmpROIFullPath);
	
	// calculate the area of each mito type (fregemented, intermidiate or elongated)
	AnalyzeMito(roiId);	
	ClearRoiManager();
	roiManager("Open", tmpROIFullPath);
}

function GenerateOverlayImages()
{
	cellsRoiPath = gImagesResultsSubFolder+"/"+gCellsRois;
	if(gManualRoi)
		cellsRoiPath += "_Manual.zip";
	else 
		cellsRoiPath += ".zip";
	ldsRoiPath = gImagesResultsSubFolder+"/"+gLDRois+".zip";
	clustersRoiPath = gImagesResultsSubFolder+"/"+gClustersRois + ".zip";
	mitoRoiPath = gImagesResultsSubFolder+"/"+gMitoRois + ".zip";
	

	//1st image: cells, LDs, and clusters ROIs in three colors on top of LDs channel
	ClearRoiManager();
	
	if(File.exists(cellsRoiPath))
	{
		roiManager("Open", cellsRoiPath);
		roiManager("Deselect");		
		roiManager("Set Color", gCellsColor);
		roiManager("Set Line Width", gRoiLineSize);
	}
	else
		print("Warning: No cells Rois found");
	nBefore = roiManager("count");
	if(File.exists(ldsRoiPath))
	{
		roiManager("Open", ldsRoiPath);
		nAfter = roiManager("count");
		for(i=nBefore;i<nAfter;i++)
		{
			roiManager("select", i);
			roiManager("Set Color", gLDsColor);
			roiManager("Set Line Width", 1);
		}
	}
	else
		print("Warning: No LDs Rois found");
		
	nBefore = roiManager("count");
	if(File.exists(clustersRoiPath))
	{
		roiManager("Open", clustersRoiPath);
		nAfter = roiManager("count");
		for(i=nBefore;i<nAfter;i++)
		{
			roiManager("select", i);
			roiManager("Set Color", gClustersColor);
			roiManager("Set Line Width", 1);
		}
	}
	else
		print("Warning: No clusters Rois found");
		
	selectImage(gLDZimageId);
	run("Enhance Contrast", "saturated=0.35");
	roiManager("Deselect");
	roiManager("Show All without labels");
	saveAs("Tiff",gImagesResultsSubFolder+"/LDs.tif");
	run("Flatten");
	saveAs("Jpeg", gImagesResultsSubFolder+"/LDs.jpg");
	
	//2nd image: cells, and mito ROIs in two colors on top of Mito. channel
	ClearRoiManager();
	
	if(File.exists(cellsRoiPath))
	{
		roiManager("Open", cellsRoiPath);
		roiManager("Deselect");		
		roiManager("Set Color", gCellsColor);
		roiManager("Set Line Width", gRoiLineSize);
	}

	nBefore = roiManager("count");
	if(File.exists(mitoRoiPath))
	{
		roiManager("Open", mitoRoiPath);
		nAfter = roiManager("count");
		for(i=nBefore;i<nAfter;i++)
		{
			roiManager("select", i);
			roiManager("Set Color", gMitoColor);
			roiManager("Set Line Width", 1);
		}
	}
	else
		print("Warning: No mito Rois found");
		
	selectImage(gMitoZimageId);
	run("Enhance Contrast", "saturated=0.35");
	roiManager("Deselect");
	roiManager("Show All without labels");
	saveAs("Tiff",gImagesResultsSubFolder+"/Mito.tif");
	run("Flatten");
	saveAs("Jpeg", gImagesResultsSubFolder+"/Mito.jpg");
}

function AnalyzeMito(roiId)
{
	//remove all Mito outside of cell
	selectImage(gMitoMaskImageId);
	//waitForUser("mito:"+getTitle());
	run("Duplicate...", "title=Cell_Mito_"+roiId+" ignore");
	roiManager("select", roiId);
	run("Clear Outside");
	//waitForUser("check mito: "+roiId);
	// turn the Mito mask into rois (2 is the value of the mito in the mask and 1 is the background)
	ClearRoiManager();	
	setThreshold(2, 1000000000000000000000000000000.0000);
	run("Analyze Particles...", "size="+0+"-Infinity display summarize add composite"); 	
	//calculate relative area of each each mito type (fregemented, intermidiate or elongated)
	roiManager("deselect");
	run("Clear Results");
	roiManager("measure");
	n = roiManager("count");
	//waitForUser("check mito n: "+n);
	totalArea = 0; fregmenetArea = 0; elongatedArea = 0; hyperElongatedArea = 0;
	for(i=0;i<n;i++)
	{
		area = Table.get("Area",i, "Results") * pixelWidth * pixelHeight;
		totalArea += area;
		if(area <  iMitoMaxFregmented)
			fregmenetArea += area;
		else if(area > iMitoMinHyperElongated)
			hyperElongatedArea += area;
		else if(area > iMitoMinElongated)
			elongatedArea += area;
	}	
	// add to cell table
	SetCompositeTables("Mito. Total area",roiId,totalArea);
	SetCompositeTables("Mito. fregmented size area (<"+iMitoMaxFregmented+"m^2)",roiId,fregmenetArea);
	SetCompositeTables("Mito. intermidiate size area",roiId,totalArea-fregmenetArea-elongatedArea-hyperElongatedArea);
	SetCompositeTables("Mito. elongated size area (>"+iMitoMinElongated+"m^2)",roiId,elongatedArea);
	SetCompositeTables("Mito. hyper elongated size area (>"+iMitoMinHyperElongated+"m^2)",roiId,hyperElongatedArea);
}
function SetCompositeTables(colName,rowId,colValue)
{
	Table.set(colName,rowId,colValue,gCompositeTable);
	if(!matches(gProcessMode, "singleFile"))
		Table.set(colName,gAllCompositeResults+rowId,colValue,gAllCompositeTable);
}

function AnalyzeLDs(roiId)
{
	//remove all LDs outside of cell
	selectImage(gLDsLabelsImageId);
	run("Duplicate...", "title=Cell_LDs_"+roiId+" ignore");	
	roiManager("select", roiId);
	run("Clear Outside");
	ClearRoiManager();
	//erode LDs labels
	erodedLabelsImageId = ErodeLabels(getImageID());
	//count the total number of LDs in cell
	run("Label image to ROIs");
	numLDs = roiManager("count");
	ClearRoiManager();
	//run("Threshold...") for SSIDC clustering
	binaryImageId = RunThreshold(1, 65535);
	// SSIDC clustring
	run("SSIDC Cluster Indicator", "distance="+iLDclusterDistance +" mindensity="+iLDminDensity);
	// remove all clusters and replace them with a single roi combining the all
	numClusters = roiManager("count");
	if(numClusters > 1)
	{
		roiManager("select", Array.getSequence(numClusters));
		roiManager("combine");
		ClearRoiManager();
		roiManager("Add");
	}
	if(numClusters > 0)
	{
		//in the labeled image of LDs in cell remove all non-clustered LDs
		selectImage(erodedLabelsImageId);
		//save clusters rois
		if(!gFirstClusterRoi)
		{
			// add prev clusters
			roiManager("open", gImagesResultsSubFolder+"/"+gClustersRois + ".zip");
		}
		else
			gFirstClusterRoi = false;

		roiManager("select", 0);
		roiManager("rename", "Cell_"+(roiId+1)+"_Clusters");
		roiManager("deselect");
		roiManager("save", gImagesResultsSubFolder+"/"+gClustersRois + ".zip");
		roiManager("select", 0);
		run("Clear Outside");
		//count the total number of clusterd LDs in cell
		ClearRoiManager()	;
		run("Label image to ROIs");
		numClusteredLDs = roiManager("count");
	}
	else {
		numClusteredLDs = 0;
	}
	
	// add to cell table
	SetCompositeTables("No. LDs",roiId,numLDs);
	SetCompositeTables("No. clusters",roiId,numClusters);
	SetCompositeTables("No. Clustered LDs",roiId,numClusteredLDs);
	
	ClearRoiManager()	;	
}

function ClearRoiManager()
{
	roiManager("reset");
}

function ErodeLabels(imageId)
{
	selectImage(imageId);
	// erode labels
	image1 = getTitle();
	Ext.CLIJ2_push(image1);
	image2 = image1+"_erode_labels";
	radius = 1.0;
	relabel_islands = false;
	Ext.CLIJ2_erodeLabels(image1, image2, radius, relabel_islands);
	Ext.CLIJ2_pull(image2);
	return getImageID();
}
function GetLDsLabeledImage(imageId)
{
	//prepare image for strardist
	gLDZimageId = DupChannelAndZProject(imageId,gLDChannel);

	//save rois and clear roi table	
	tmpROIFullPath = gImagesResultsSubFolder+"/"+gTempROIFile;
	roiNotEmpty = SaveROIs(tmpROIFullPath);
	//after saving rois clear roi table
	ClearRoiManager();	
	
	GetLDsRois();
	
	ClearRoiManager();
	
	if(roiNotEmpty)
		roiManager("Open", tmpROIFullPath);
	
	return gLDsLabelsImageId;	
}

function ScaleImage(imageId, scaleFactor)
{
	s = "x="+scaleFactor
	+" y="+scaleFactor
	//+" width="+width*scaleFactor
	//+" height="+height*scaleFactor
	+" interpolation=None create";
	run("Scale...",s);
	return getImageID();
}

function SaveROIs(fullPath)
{
	if(roiManager("count") <= 0)
		return false;
	roiManager("deselect");
	roiManager("save", fullPath);	
	return true;
}

function UsePrevRun(title,usePrevRun,restoreRois)
{
	if(!usePrevRun)
		return false;

	savePath = gImagesResultsSubFolder + "/"+gSaveRunDir+"/";
	labeledImageFullPath = savePath + title +".tif";
	labeledImageRoisFullPath = savePath + title +"_RoiSet.zip";

	if(!File.exists(labeledImageFullPath))
		return false;
	
	if(restoreRois && !File.exists(labeledImageRoisFullPath))
		return false;
		
	open(labeledImageFullPath);
	rename(title);
	id = getImageID();
	if(restoreRois)
		openROIs(labeledImageRoisFullPath,true);
	selectImage(id);
	//print("Using stored "+title+" labeled image");
	return true;;
}


function StoreRun(title,storeROIs)
{
	//waitForUser("store title: "+title);
	savePath = gImagesResultsSubFolder + "/"+gSaveRunDir+"/";
	//waitForUser("store gImagesResultsSubFolder: "+gImagesResultsSubFolder);
	//waitForUser("store savePath: "+savePath);
	
	labeledImageFullPath = savePath + title +".tif";
	
	selectWindow(title);
	saveAs("Tiff", labeledImageFullPath);
	rename(title);
	if(storeROIs)
	{
		labeledImageRoisFullPath = savePath + title +"_RoiSet.zip";
		SaveROIs(labeledImageRoisFullPath);
	}
}

function RunStarDistModel(imageId)
{
	labelImageId = -1;


//	labeledImageFullPath = gImagesResultsSubFolder + "/" + StarDistWindowTitle +".tif";
//	labeledImageRoisFullPath = gImagesResultsSubFolder + "/" + StarDistWindowTitle +"_RoiSet.zip";

	if(UsePrevRun(gStarDistWindowTitle,iUseStarDistPrevRun,true))
	{
		print("Using StarDist stored labeled image");
		labelImageId = getImageID();
	}
	else 
	{
		starDistModel = "'Versatile (fluorescent nuclei)'";
		nmsThresh = 0.4;
		
		print("Progress Report: StarDist started. That might take a few minutes");	
		selectImage(imageId);
		title = getTitle();
		run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=["
		+"'input':'"+title+"'"
		+", 'modelChoice':"+starDistModel
		+", 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.8'"
		+", 'probThresh':'"+iSDprobThreshold+"'"
		+", 'nmsThresh':'"+nmsThresh+"'"
		+", 'outputType':'Both', 'nTiles':'4', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");		

		labelImageId = getImageID();
		StoreRun(gStarDistWindowTitle,true);
		print("Progress Report: StarDist ended.");	
		if(roiManager("count") <= 0)
		{
			print("WARNING!!!: Stardist did not identify any cell/object in " + title);
		}
		print("num lds: "+roiManager("count"));
	}
	return labelImageId;
}
function RunIlastikModel(imageId)
{
	selectImage(imageId);
	title = getTitle();
	found = false;
	IlastikSegmentationOutFile = title+gIlastikSegmentationExtention;
	IlastikOutFilePath = gImagesResultsSubFolder+"/"+gSaveRunDir+"/";
	if (iUseIlastikPrevRun)
	{
		if (File.exists(IlastikOutFilePath+IlastikSegmentationOutFile))
		{
			print("Reading existing Ilastik AutoContext output ...");
			//run("Import HDF5", "select=[A:\yairbe\Ilastic Training\Cre off HD R.h5] datasetname=/data axisorder=tzyxc");
			//run("Import HDF5", "select=["+resFolderSub+IlastikSegmentationOutFile+"] datasetname=/exported_data axisorder=yxc");
			run("Import HDF5", "select=["+IlastikOutFilePath+IlastikSegmentationOutFile+"] datasetname=/data axisorder=tzyxc");

			//rename("Segmentation");
			rename(IlastikSegmentationOutFile);
						
			found = true;
		}
	}
	if (!found)
	{
		print("Progress Report: Ilastik AutoContext classifier started. That might take a few minutes");	
		//run("Run Autocontext Prediction", "projectfilename=[A:\\yairbe\\Ilastic Training\\CreOFF-Axon-Classifier_v133post3.ilp] 
		//    inputimage=[A:\\yairbe\\Ilastic Training\\Cre off HD R.h5\\data] autocontextpredictiontype=Segmentation");

		run("Run Autocontext Prediction", "projectfilename=["+iIlastikModelPath+"] inputimage=["+title+"] autocontextpredictiontype=Segmentation");
		//rename("Segmentation");
		rename(IlastikSegmentationOutFile);

		// save Ilastik Output File
		selectWindow(IlastikSegmentationOutFile);
		print("Saving Ilastik autocontext classifier output...");
		//run("Export HDF5", "select=["+resFolder+IlastikProbOutFile1+"] exportpath=["+resFolder+IlastikProbOutFile1+"] datasetname=data compressionlevel=0 input=["+IlastikProbOutFile1+"]");	
		run("Export HDF5", "select=["+IlastikOutFilePath+IlastikSegmentationOutFile+"] exportpath=["+IlastikOutFilePath+IlastikSegmentationOutFile+"] datasetname=data compressionlevel=0 input=["+IlastikSegmentationOutFile+"]");	
		print("Progress Report: Ilastik ended.");	
	}	
	rename(IlastikSegmentationOutFile);
	//setVoxelSize(width, height, depth, unit); multiplying area size instead
	return getImageID();
}

function RunCellposeModel(imageId, cellposeModel)
{
	isDefaultModel = false;
	if(cellposeModel == gCadmiumCellposeModel)
	{
		isDefaultModel = true; cellposeModelPath = "cyto2"; useCellposePrevRun = iUseCellposePrevRun; cellposeCellDiameter = iParticleCellposeCellDiameter; cellposeProbThreshold = 0.0; cellposeFlowThreshold = 0.4;
	}
	else
	{
		print("Error: Unidentified Cellpose model: "+cellposeModel);
		return -1;
	}
	selectImage(imageId);
	title = getTitle();
	CellposeWindowTitle = "label image - ";
	
	//if cellpose was used externaly to generate the label map of the cells
	//it will be stored in the input directory under Cellpose/Segmentation
	if(UseExternalRun("Cellpose", CellposeWindowTitle, cellposeModel))
	{
		print("Using "+cellposeModel+ "Cellpose external run generated labeled image");
		labelImageId = getImageID();
	}
	else if(UsePrevRun(CellposeWindowTitle,useCellposePrevRun,false))
	{
		print("Using "+cellposeModel+ " Cellpose stored labeled image");
		labelImageId = getImageID();
	}
	else 
	{
		print("Progress Report: "+cellposeModel+ " started. That might take a few minutes");	
		if(isDefaultModel)
		{
			parms = "diameter="+cellposeCellDiameter
					+" cellproba_threshold="+cellposeProbThreshold
					+" flow_threshold="+cellposeFlowThreshold
					+" anisotropy=1.0 diam_threshold=12.0"
					+" model="+cellposeModelPath
					+" nuclei_channel=0 cyto_channel=0 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=";
			print("parms:"+parms);
			print("end of parms");
			run("Cellpose Advanced", parms);
		}
		else
			run("Cellpose Advanced (custom model)", "diameter="+cellposeCellDiameter
				+" cellproba_threshold="+cellposeProbThreshold
				+" flow_threshold="+cellposeFlowThreshold
				+" anisotropy=1.0 diam_threshold=12.0"
				+" model_path="+File.getDirectory(cellposeModel)
				+" model="+cellposeModelPath
				+" nuclei_channel=0 cyto_channel=0 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additionnal_flags=");
		labelImageId = getImageID();
		rename(CellposeWindowTitle);
		//waitForUser("title:"+title);
		StoreRun(CellposeWindowTitle,false);
		print("Progress Report: "+cellposeModel+ " ended.");	
	}	

	return labelImageId;
}

function UseExternalRun(app, title, cellposeModel)
{
	if(app == "Cellpose")
	{
		subDir = gCellposeExtRunSubdir;
		fileSuffix = gCellposeExtRunFileSuffix;
	}
	else{
		print("Warning: Wrong app: " + app + ". Ignored");
		return false;
	}
	labeledImageFullPath = File.getDirectory(gFileFullPath)+cellposeModel+subDir+gFileNameNoExt+fileSuffix;

	if(File.exists(labeledImageFullPath))
	{
		open(labeledImageFullPath);
		rename(title);
		id = getImageID();
		selectImage(id);	
		return true;	
	}
	return false;
}

function GenerateROIsFromLabelImage(labelImageId,type,filterMinAreaSize)
{
	nBefore = roiManager("count");
	selectImage(labelImageId);
	run("Label image to ROIs");
	nAfter = roiManager("count");
	if(filterMinAreaSize > 0)
	{
		run("Clear Results");
		roiManager("measure");
		for(i=nAfter-1;i>=nBefore;i--)
		{
			area = Table.get("Area",i, "Results");
			if(area < filterMinAreaSize)
			{
				roiManager("select", i);
				roiManager("delete");
				nAfter--;
			}
		}
	}
	if(type != "")
	{
		for(i=nBefore;i<nAfter;i++)
		{
			roiManager("select", i);
			roiManager("rename", type+"_"+(i-nBefore+1));
		}
	}
	return nAfter - nBefore;
}
function DupChannelAndZProject(imageId, channel)
{
	selectImage(imageId);
	//1. prepare image for Cellpose: duplicate the 3rd channel and z-project it
	run("Duplicate...", "duplicate channels="+channel);
	run("Z Project...", "projection=[Max Intensity]");
	return getImageID();	
}

function RemoveArtifacts(imgId,spotThreshold,channel)
{
	roiManager("Deselect");
	run("Clear Results");
	run("Select None");
	
	selectImage(imgId);	
	Stack.setChannel(channel);
	run("Set Measurements...", "mean min display redirect=None decimal=3");
	run("Measure");
	meanIntensity = Table.get("Mean",0, "Results");
	run("Macro...", "code=if(v>" + spotThreshold + ")v=" + meanIntensity); //Notice: no spaces allowed in macro
	run("Clear Results");
}

function StoreROIs(path,fileName)
{
	SaveROIs(path +"/" + fileName+".zip");
}
function compositeAreas(imageId)
{
	selectImage(imageId);

	// set the B/C of CD35 and CD25 channels
	Stack.setChannel(CD35_CHANNEL);
	setMinAndMax(0, iCD35MaxPixelDisplayValue);
	Stack.setChannel(CD23_CHANNEL);
	setMinAndMax(0, iCD23MaxPixelDisplayValue);
	
	//make a composite image of the two channels
	channels = newArray(4);
	channels[CD35_CHANNEL-1] = 1;
	channels[CD23_CHANNEL-1] = 1;	
	Stack.setDisplayMode("composite");
	Stack.setActiveChannels(String.join(channels,""));

	
	run("Select None");
	roiManager("Show None");
	
	// show CD23, CD35 and GC on the composite image and flatten them on it
	rois = newArray(SelectRoiByName("GC"),SelectRoiByName("CD35"),SelectRoiByName("CD23"));
	colors = newArray(iGCColor,iCD35Color,iCD23Color);

	prevImgId = 0;
	for(i=0;i<rois.length;i++)
	{
		roiManager("Select", rois[i]);
		RoiManager.setGroup(0);
		RoiManager.setPosition(0);	
		roiManager("Set Color", colors[i]);
		roiManager("Set Line Width", iROILineWidth);
		prevImgId = imageId;
		run("Flatten");	
		imageId = getImageID();
		if(i > 0)
		{
			selectImage(prevImgId);
			close();
			selectImage(imageId);						
		}
	}
	saveAs("Jpeg", gImagesResultsSubFolder+"/"+"CD35nCD23."+"jpeg");
	return imageId;
}

function StroeTCellsInfo(CD35Area, CD23Area)
{
	GenerateLineHistogram(gTCellsBitmapImgId,"T Cells Count", false,true);
	tableName = "T Cells info";
	TCellsInCD35 = CountWhitePixelsInROI(gTCellsBitmapImgId, SelectRoiByName("CD35"));
	TCellsInCD23 = CountWhitePixelsInROI(gTCellsBitmapImgId, SelectRoiByName("CD23"));
	Table.create(tableName);
	Table.set("CD35 Count", 0,TCellsInCD35,tableName);
	Table.set("CD35 Density", 0,d2s(TCellsInCD35/CD35Area,5),tableName);
	Table.set("CD23 Count", 0,TCellsInCD23,tableName);
	Table.set("CD23 Density", 0,d2s(TCellsInCD23/CD23Area,5),tableName);
	Table.set("CD35 Basal Count", 0,(TCellsInCD35 - TCellsInCD23),tableName);
	Table.set("CD35 Basal Density", 0,d2s((TCellsInCD35 - TCellsInCD23)/(CD35Area- CD23Area),5),tableName);
  	Table.save(gImagesResultsSubFolder + "/"+tableName+".csv", tableName);
}

function GenerateLineHistogram(imageId,title, normalizeY, convertToMicrons)
{

	selectImage(imageId);
	SelectRoiByName("Line");
	run("Plot Profile");	
	saveAs("Jpeg", gImagesResultsSubFolder+"/"+title+"_LineHistogram.jpg");
	//now save as csv file as well
  	Plot.getValues(x, y);
  	if(convertToMicrons)
  	{
		for(i=0;i<x.length;i++)
			x[i] *= pixelWidth;
  	}
  	run("Close") ;
  	tableName = title+"_Histogram_values";
  	Array.show(tableName, x, y);
  	//Table.save(gResultsSubFolder+"/"+tableName+".csv", tableName);
  	Table.save(gImagesResultsSubFolder + "/"+tableName+".csv", tableName);
  	run("Close") ;
  	//save normalized table as well
  	normalizeXArray(x,iLineMargin);
  	if(normalizeY)
 		normealizeYArray(y,title);
   	//Table.save(gResultsSubFolder+"/"+tableName+".csv", tableName);
  	tableName = title+"_Histogram_values_normalized";
  	Array.show(tableName, x, y);
  	Table.save(gImagesResultsSubFolder + "/"+tableName+".csv", tableName);
	run("Close") ;
}


function CountWhitePixelsInROI(BitmapImgId,RoiId)
{
	selectImage(BitmapImgId);
	run("Set Measurements...", "area limit display redirect=None decimal=3");
	roiManager("Select", RoiId);
	setAutoThreshold("Default dark no-reset");
	setThreshold(1, 255);
	run("Measure");
	run("Clear Results");
	roiManager("Measure");
	count = Table.get("Area",0, "Results");
	return count;
}


//FlattenRois:
//to a given image it flattens all rois onto it and stroes it as Jpeg
function FlattenRois(imageId, name, color, width, extention,withLabels)
{
	selectImage(imageId);
	
	//RoiManager.setPosition(3);
	roiManager("Deselect");
	roiManager("Set Color", color);
	roiManager("Set Line Width", width);
	//n = roiManager("count");
	//arr = Array.getSequence(n);
	//roiManager("select", arr);
	//Roi.setFontSize(5)
	if(withLabels)
		roiManager("Show All with labels");
	else 
		roiManager("Show All without labels");
	//	Overlay.setLabelFontSize(5, 'sacle')
	run("Flatten");
	saveAs(extention, gImagesResultsSubFolder+"/"+name+"_rois."+extention);
}

// CreateSegmentationAndDensityMap:
// 1. on a given channel, run Tubeness filter + Threshold to generate Mask
// 2. run mean intensisty on the result to get density map
// 3. store it to disk 
// 4. stroe to disk the original channel with the threshold mask
function CreateSegmentationAndDensityMap(window, minThreshold, minAxonSize)
{
	selectWindow(window);
	getPixelSize(unit,pixelWidth, pixelHeight);
	pixelDensityRadius = Math.round((iDensityRadius)/pixelWidth);
	//print("pixelDensityRadius: " + pixelDensityRadius);
	run("Tubeness", "sigma="+iTubenessSigma+" use");
	setThreshold(minThreshold, 1000000000000000000000000000000.0000);
	setOption("BlackBackground", false);
	run("Analyze Particles...", "size="+minAxonSize+"-Infinity show=Masks exclude composite");
	run("Convert to Mask");	
	run("Create Selection");
	run("Select None");
	
	run("Divide...", "value=2.55");
	//run("Brightness/Contrast...");
	//setMinAndMax(0, iMaxDensity);

	run("Mean...", "radius="+pixelDensityRadius);
	run("Fire");
	setMinAndMax(0, iMaxDensity);	
	run("Calibration Bar...", "location=[Upper Right] fill=White label=Black number=5 decimal=0 font=12 zoom="+iZoomLevel+" overlay");
	dansityMapTitle = gFileNameNoExt + "_" + window + "_DensityMap_R"+pixelDensityRadius;
	rename(dansityMapTitle);
	//saveAs("Tiff", gImagesResultsSubFolder+"/"+title+".tif");
	saveAs("Tiff", gImagesResultsSubFolder+"/"+window+".tif");
	dansityMapTitle = getTitle();
	
	selectWindow(window);
	run("Grays");
	run("Restore Selection");
	//saveAs("Tiff", gImagesResultsSubFolder+"/"+gFileNameNoExt + "_"+window+"_Segmentation_T"+iMinThreshold+".tif");
	saveAs("Tiff", gImagesResultsSubFolder+"/"+window+"_Segmentation_T"+minThreshold+".tif");
	return dansityMapTitle;
}

function FinalActions()
{
	if(gAllCompositeResults > 0) // stroe allCompositeTable table
		Table.save(gResultsSubFolder+"/"+gAllCompositeTable+".csv", gAllCompositeTable);
}
// end of single file analysis

//--------Helper functions-------------

function Initialization()
{
	requires("1.53c");
	run("Check Required Update Sites");
	// for CLIJ
	run("CLIJ2 Macro Extensions", "cl_device=");
	Ext.CLIJ2_clear();

	//run("Configure ilastik executable location", "executablefile=["+iIlastikExe+"] numthreads=-1 maxrammb=150000");
	run("Cellpose setup...", "cellposeenvdirectory="+iCellposeEnv+" envtype=conda usegpu=true usemxnet=false usefastmode=false useresample=false version=2.0");		
	
	setBatchMode(false);
	run("Close All");
	close("\\Others");
	print("\\Clear");
	run("Options...", "iterations=1 count=1 black");
	run("Set Measurements...", "area display redirect=None decimal=3");
	roiManager("Reset");

	CloseTable("Results");
	CloseTable(gCompositeTable);	
	CloseTable(gAllCompositeTable);

	run("Collect Garbage");

	if (gBatchModeFlag)
	{
		print("Working in Batch Mode, processing without opening images");
		setBatchMode(gBatchModeFlag);
	}	

}

function checkInput()
{
	getDimensions (ImageWidth, ImageHeight, ImageChannels, ImageSlices, ImageFrames);

	if(ImageChannels < 3)
	{
		print("Fatal error: input file must include 3 channels: Lipid Droplets, Litosomes, and Mitocondria stainings");
		return false;
	}
	getPixelSize(unit,pixelWidth, pixelHeight);
	if(!matches(unit, "microns") && !matches(unit, "um"))
	{
		print("Fatal error. File " + gFileFullPath + " units are "+ unit+ " and not microns");
		return false;
	}
	return true;
}
//------openROIsFile----------
//open ROI file with 
function openROIsFile(ROIsFileNameNoExt, clearROIs)
{
	roiManager("deselect");
	// first delete all ROIs from ROI manager
	if(clearROIs && roiManager("count")
		roiManager("delete");

	// ROIs are stored in "roi" suffix in case of a single roi and in "zip" suffix in case of multiple ROIs
	RoiFileName = ROIsFileNameNoExt+".roi";
	ZipRoiFileName = ROIsFileNameNoExt+".zip";
	if (File.exists(RoiFileName) && File.exists(ZipRoiFileName))
	{
		if(File.dateLastModified(RoiFileName) > File.dateLastModified(ZipRoiFileName))
			roiManager("Open", RoiFileName);
		else
			roiManager("Open", ZipRoiFileName);
		return true;
	}
	if (File.exists(RoiFileName))
	{
		roiManager("Open", RoiFileName);
		return true;
	}
	if (File.exists(ZipRoiFileName))
	{
		roiManager("Open", ZipRoiFileName);
		return true;
	}
	return false;
}

function openROIs(ROIsFullName, clearROIs)
{
	roiManager("deselect");
	// first delete all ROIs from ROI manager
	if(clearROIs && roiManager("count") > 0)
		roiManager("delete");

	if (File.exists(ROIsFullName))
	{
		roiManager("Open", ROIsFullName);
		return true;
	}
	return false;
}



//----------LoopFiles-------------
// according to gProcessMode analyzes a single file, or loops over a directory or sub-directories
function LoopFiles()
{
	//SetProcessMode();
	gMainDirectory = File.getParent(iCellposFile);
	gFileFullPath = gMainDirectory;
	gResultsSubFolder = gMainDirectory + File.separator + iResultsDir + File.separator; 
	File.makeDirectory(gResultsSubFolder);
	SaveParms(gResultsSubFolder);
	return ProcessFile(gMainDirectory); 	
	return true;
}
function SaveParms(resFolder)
{
	//waitForUser("macro"+File.getNameWithoutExtension(getInfo("macro.filepath")));
	// print parameters to Prm file for documentation
	PrmFile = pMacroName+"Parameters.txt";
	if(gProcessMode == "singleFile")
		PrmFile = resFolder + File.getNameWithoutExtension(gFileFullPath) + "_" + PrmFile;
	else 
		PrmFile = resFolder + PrmFile;
		
	File.saveString("macroVersion="+pMacroVersion, PrmFile);
	File.append("", PrmFile); 
	
	File.append("RunTime="+getTimeString(), PrmFile);
	
	// save user input
	File.append("processMode="+gProcessMode, PrmFile);  
	File.append("iCellposFile="+iCellposFile+" \n", PrmFile); 
	File.append("iParticlesFile="+iParticlesFile+" \n", PrmFile); 
	File.append("iMnFile="+iMnFile, PrmFile)
	File.append("iMnMinThreshold="+iMnMinThreshold, PrmFile)
 	//File.append("iparticleMinThreshold="+iparticleMinThreshold, PrmFile)
 	File.append("iResultsDir="+iResultsDir, PrmFile)
	File.append("iCellposeEnv="+iCellposeEnv, PrmFile)
	File.append("iParticleCellposeCellDiameter="+iParticleCellposeCellDiameter, PrmFile)
	File.append("iUseCellposePrevRun="+iUseCellposePrevRun, PrmFile)

}
function getTimeString()
{
	MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
	DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	TimeString ="Date: "+DayNames[dayOfWeek]+" ";
	if (dayOfMonth<10) {TimeString = TimeString+"0";}
	TimeString = TimeString+dayOfMonth+"-"+MonthNames[month]+"-"+year+", Time: ";
	if (hour<10) {TimeString = TimeString+"0";}
	TimeString = TimeString+hour+":";
	if (minute<10) {TimeString = TimeString+"0";}
	TimeString = TimeString+minute+":";
	if (second<10) {TimeString = TimeString+"0";}
	TimeString = TimeString+second;
	return TimeString;
}
//===============================================================================================================
// Loop on all files in the folder and Run analysis on each of them
function ProcessFiles(directory) 
{
	Table.create(gAllCompositeTable);		
	gAllCompositeResults = 0;

	setBatchMode(gBatchModeFlag);
	dir1=substring(directory, 0,lengthOf(directory)-1);
	idx=lastIndexOf(dir1,File.separator);
	subdir=substring(dir1, idx+1,lengthOf(dir1));

	// Get the files in the folder 
	fileListArray = getFileList(directory);
	
	// Loop over files
	for (fileIndex = 0; fileIndex < lengthOf(fileListArray); fileIndex++) {
		if (endsWith(fileListArray[fileIndex], iFileExtension) ) {
			gFileFullPath = directory+File.separator+fileListArray[fileIndex];
			print("\nProcessing:",fileListArray[fileIndex]);
			showProgress(fileIndex/lengthOf(fileListArray));
			if(!ProcessFile(directory))
				return false;
			CleanUp(false);		
		} // end of if 
	} // end of for loop
	FinalActions();
	CleanUp(true);
	return true;
} // end of ProcessFiles

function CleanUp(finalCleanUp)
{
	run("Close All");
	close("\\Others");
	run("Collect Garbage");
	if (finalCleanUp) 
	{
		CloseTable(gAllCompositeTable);	
		setBatchMode(false);
	}
}
function SetProcessMode()
{
		// Choose image file or folder
	if (matches(gProcessMode, "singleFile")) {
		gFileFullPath=File.openDialog("Please select an image file to analyze");
		gMainDirectory = File.getParent(gFileFullPath);
	}
	else if (matches(gProcessMode, "wholeFolder")) {
		gMainDirectory = getDirectory("Please select a folder of images to analyze"); }
	
	else if (matches(gProcessMode, "AllSubFolders")) {
		gMainDirectory = getDirectory("Please select a Parent Folder of subfolders to analyze"); }
}

//===============================================================================================================
function CloseTable(TableName)
{
	if (isOpen(TableName))
	{
		selectWindow(TableName);
		run("Close");
	}
}
