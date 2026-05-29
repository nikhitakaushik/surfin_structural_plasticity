import pickle 
import pandas as pd
import numpy as np

def get_single_day_df(fname,mouse_id,fov_id):
    """
    Load single day pickle file & return spine-level DataFrame. 

    Extracts spine positions, volumes, flags, & groupings from a pickle data object. 
    Assigns unique spine & dendrite IDs.

    Parameters
    ----------
    fname : str
        Full path to pickle file. 
    mouse_id : str
        Mouse identifier (ex. JL117) used to name spine/dendrite IDs
    fov_id : str
        FOV identifier (ex. FOV1) used to name spine/dendrite IDs/ 

    Returns
    -------
    df : pd.DataFrame
        Columns: spine_ID, dendrite_ID, position, volume, flags. 
    """
    with open(fname,"rb") as f:
        data = pickle.load(f)
    
    #DELETE
    #print("spine_groupings:", data.parameters["Spine Groupings"])
    #print("num positions:", len(data.ROI_positions["Spine"]))

    if data.parameters['Spine Groupings'] == []:
        spine_groupings = [(list(range(len(data.ROI_positions["Spine"]))))]
    else:
        spine_groupings = data.parameters["Spine Groupings"]
    spine_positions = data.ROI_positions["Spine"]
    spine_volumes = data.corrected_spine_volume
    spine_flags = data.ROI_flags["Spine"]

    spine_params = {"spine_ID": [], "dendrite_ID": []}

    dendrites = {}
    for i, dendrite in enumerate(spine_groupings):
        dendrites[f'dendrite_{i+1}'] = dendrite

    for index,x in enumerate(spine_positions):
        if any(index in val for val in dendrites.values()):
            dendrite_spine = next(key for key, val in dendrites.items() if index in val)
            index = str(index)
            spine_id = mouse_id + "_" + fov_id + "_" + dendrite_spine + "_" + index
            dendrite_spine = mouse_id + "_" + fov_id + "_" + dendrite_spine
            spine_params["spine_ID"].append(spine_id)
            spine_params["dendrite_ID"].append(dendrite_spine)

    spine_params["position"] = spine_positions
    spine_params["volume"] = spine_volumes
    spine_params["flags"] = spine_flags

    # DELETE
    # for key, val in spine_params.items():
        # print(f"  {key}: {len(val) if hasattr(val, '__len__') else 'scalar'}")
    df = pd.DataFrame(spine_params)

    return df

def find_max_pos_diff(x):
    """
    Calculate maximum pairwise distance between spine positions. 

    Used to estimate dendrite length from the furthest apart pair of spine positions along a single dendrite. 

    Parameters
    ----------
    x : np.ndarray
        1D array of spine positions along a dendrite. 

    Returns
    -------
    max : float
        Maximum pairwise distance between any two positions in x. 
    """
    pair_dif = x[:, np.newaxis] - x[np.newaxis, :]
    max = np.nanmax(pair_dif)
    return max
    
