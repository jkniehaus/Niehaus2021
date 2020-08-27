macro "autoAdjust" {
        /*
         * rewriting "Auto contrast adjustment" button of "Brightness/Contrast"
         * Damien Guimond
         * 20120516
         * Acknowledgements: Kota Miura
         */
       
        AUTO_THRESHOLD = 5000;
       
        getRawStatistics(count, mean, min, max, std);
        limit = count/10;
        threshold = count/AUTO_THRESHOLD;
       
        nBins = 256;
        getHistogram(values, counts, nBins);
       
        i = -1;
        condition = 0;
        while (condition != 1) {
                i = i+1;
                count = counts[i];
                if (count > limit) { count = 0; }
                found = count > threshold;
                if (found) { condition = 1; }
                else if (i >= 255) { condition = 1; }
        }
        hmin = i;
       
        i = 256;
        condition = 0;
        while (condition != 1) {
                i = i-1;
                count = counts[i];
                if (count > limit) { count = 0; }
                found = count > threshold;
                if (found) { condition = 1; }
                else if (i < 1) { condition = 1; }
        }
        hmax = i;
       
        setMinAndMax(hmin, hmax);
        run("Apply LUT");
}
