//Start by installing the autoAdjust macro and adding an outline of the ipsilateral DH (or other region of interest) to ROI manager
run("Select None");
run("Duplicate...", "duplicate channels=2");
roiManager("Select", 0);
run("Clear Outside");
run("Select None");
run("autoAdjust")

setAutoThreshold("Default dark");
run("Convert to Mask");
//analyze area
roiManager("Select", 0);
run("Analyze Particles...", "size=0-Infinity summarize");
close();
//iba1+ cells; create ROI for cells
run("Duplicate...", "duplicate channels=1");
run("autoAdjust")
roiManager("Select", 0);
run("Clear Outside");
run("Select None");
setAutoThreshold("Huang dark");
run("Convert to Mask");
run("Despeckle");
run("Watershed");
run("Create Selection");
roiManager("Add");
close();
//Count number of cells that contain x area of IBA1
run("Duplicate...", "duplicate channels=2");
roiManager("Select", 0);
run("Clear Outside");
run("Select None");
run("autoAdjust")

setAutoThreshold("Default dark");

run("Convert to Mask");
roiManager("Select", 1);
run("Clear Outside");
run("Select None");
roiManager("Select",0);
run("Analyze Particles...", "size=20-Infinity show=Nothing summarize");
roiManager("Select", 1);
roiManager("Delete");
roiManager("Select", 0);
//roiManager("Delete");
//close();
close();

