// === USER INPUT ===
mouseID = getString("Enter Mouse ID (e.g., Mouse01):", "");

// === INIT PATHS ===
basePath = "Z:/People/JenLi/Imaging_Data/L23_L5_structural_imaging/motion_corrected/";
mousePath = basePath + mouseID + "/";
logMessages = newArray();
logCount = 0;

// === FIND VALID DATE FOLDERS ===
allFolders = getFileList(mousePath);
validDates = newArray();
validDateCount = 0;

for (i = 0; i < allFolders.length; i++) {
    folder = trim(allFolders[i]);
    if (endsWith(folder, "/") && lengthOf(folder) == 7 && matches(folder, "^[0-9]{6}/$")) {
        validDates[validDateCount++] = substring(folder, 0, 6); // remove trailing slash
    }
}

if (validDateCount == 0) {
    showMessage("No valid date folders (YYMMDD) found for " + mouseID);
    exit();
}

// === FIND FOVs PRESENT IN ALL DATES ===
firstDateFOVs = getFileList(mousePath + validDates[0] + "/");
commonFOVs = newArray();

for (i = 0; i < firstDateFOVs.length; i++) {
    folder = firstDateFOVs[i];
    if (endsWith(folder, "/")) {
        fov = substring(folder, 0, lengthOf(folder) - 1);
        commonFOVs[lengthOf(commonFOVs)] = fov;
    }
}

filteredFOVs = newArray();
for (i = 0; i < lengthOf(commonFOVs); i++) {
    fov = commonFOVs[i];
    presentInAll = true;
    for (j = 1; j < validDateCount; j++) {
        fovPath = mousePath + validDates[j] + "/" + fov + "/";
        if (!File.isDirectory(fovPath)) {
            presentInAll = false;
            break;
        }
    }
    if (presentInAll) {
        filteredFOVs[lengthOf(filteredFOVs)] = fov;
    }
}
commonFOVs = filteredFOVs;

if (lengthOf(commonFOVs) == 0) {
    showMessage("No FOVs found that are present across all dates for " + mouseID);
    exit();
}

// === PROMPT FOR FOV SELECTION ===
fovListRaw = joinFOVs(commonFOVs);
selectedFOVsCSV = getString("Enter FOVs to process (comma-separated):", fovListRaw);
selectedFOVs = split(selectedFOVsCSV, ",");
for (i = 0; i < selectedFOVs.length; i++) {
    selectedFOVs[i] = trim(selectedFOVs[i]);
}

// === STACKING OPERATION ===
sliceIndices = newArray();
channelNames = newArray();
channelCount = 0;

for (f = 0; f < selectedFOVs.length; f++) {
    fovName = selectedFOVs[f];

    for (j = 0; j < validDateCount; j++) {
        date = validDates[j];
        fovPath = mousePath + date + "/" + fovName + "/";
        subdirs = getFileList(fovPath);

        if (subdirs.length == 0) {
            msg = "No channel folders in: " + fovPath;
            print(msg);
            logMessages[logCount++] = msg;
            continue;
        }

        for (s = 0; s < subdirs.length; s++) {
            channel = subdirs[s];
            if (!endsWith(channel, "/")) continue;

            tifPath = fovPath + channel + "XYZAligned.tif";
            if (!File.exists(tifPath)) {
                msg = "Missing: " + tifPath;
                print(msg);
                logMessages[logCount++] = msg;
                continue;
            }

            open(tifPath);
            origTitle = getTitle();
            run("Z Project...", "projection=[Average Intensity]");
            projTitle = "AVG_" + origTitle;
            label = date + "_" + fovName;

            channelName = substring(channel, 0, lengthOf(channel) - 1);
            stackTitle = "Stack_" + fovName + "_" + channelName;

            index = -1;
            for (k = 0; k < channelCount; k++) {
                if (channelNames[k] == stackTitle) {
                    index = k;
                    break;
                }
            }

            if (index == -1) {
                selectWindow(projTitle);
                rename(stackTitle);
                channelNames[channelCount] = stackTitle;
                sliceIndices[channelCount] = 2;

                selectWindow(stackTitle);
                setSlice(1);
                setMetadata("Label", label);

                channelCount++;
            } else {
                selectWindow(projTitle);
                run("Copy");
                selectWindow(stackTitle);
                run("Add Slice");
                run("Paste");
                setSlice(sliceIndices[index]);
                setMetadata("Label", label);
                sliceIndices[index]++;
                close(projTitle);
            }

            selectWindow(origTitle);
            close();
        }
    }
}

// === OPTIONAL: SAVE LOG FILE ===
/*
logFile = File.saveDialog("Save log file as", mouseID + "_log.txt");
if (logFile != null) {
    f = File.open(logFile);
    for (i = 0; i < logCount; i++) {
        print(f, logMessages[i]);
    }
    File.close(f);
}
*/

// === HELPER: Join FOVs into a string ===
function joinFOVs(arr) {
    out = "";
    for (i = 0; i < arr.length; i++) {
        out += arr[i];
        if (i < arr.length - 1) out += ", ";
    }
    return out;
}
