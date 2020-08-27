//Macro is run for each image; GFAP is channel 1, DAPI is channel 2
//Prior to running, install autoAdjust macro and add DH outline to ROI manager for a given image
run("Duplicate...", "duplicate channels=1");
roiManager("Select", 0);
run("Clear Outside")
run("Select None")
run("Despeckle");
run("autoAdjust")
setAutoThreshold("Default dark no-reset");
//setThreshold(78, 255);
run("Convert to Mask");
roiManager("Select",0);
run("Analyze Particles...", "summarize");
run("Select None")
close();
roiManager("Delete");
