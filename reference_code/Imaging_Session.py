import os
import json
import pickle
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

import dend_imaging_utils as d_utils
import dend_imaging_load_align as load_align
import dend_imaging_analysis as analysis
import dend_imaging_plots as plots

class ImagingSession:
    """Encapsulates one imaging experiment (single mouse + date)."""

    def __init__(self, exp_info,file_paths, reanalyze=False):
        """
        Initialize and (optionally) reprocess the session.

        Parameters
        ----------
        mouse_id : str
        date : str
        base_dir : str
            Root folder for analyzed data (e.g., "C:/Users/Jenny/Desktop/Analyzed_data/individual")
        ws_data_dir : str
            Root folder for raw wavesurfer files.
        reanalyze : bool
            Whether to re-run extraction/processing instead of loading existing files.
        """

        self.reanalyze = reanalyze

        # --- load paths
        self.img_dir = file_paths["data_dir"]    
        self.lever_dir = os.path.join(file_paths["ws_dir"], exp_info["date"], exp_info["mouse_id"])

        # # --- save paths
        # self.savepath = os.path.join(file_paths["base_dir"], exp_info["mouse_id"], exp_info["date"],
        #                             f'FOV{exp_info["fov_id"]}',"processed_data",f'cell{exp_info["cell_id"]}')

        # if not os.path.exists(self.savepath):
        #     os.makedirs(self.savepath)
        # if not os.path.exists(os.path.join(self.savepath,"analysis")):
        #     os.mkdir(os.path.join(self.savepath,"analysis"))
        #     os.mkdir(os.path.join(self.savepath,"figures"))

        # self.session_data = None
        # self.group_data = {}
        # self.results = {}
        # self.meta = None
        # self.exp_info = exp_info
        

    # -----------------------------------------------------------
    # Loading / alignment
    # -----------------------------------------------------------

    def load_align_dff_lever(self,reanalyze=False):
        """Load existing data or perform full extraction."""
        data_path = os.path.join(self.savepath, "session_data.pkl")
        meta_path = os.path.join(self.savepath, "metadata.json")

        if reanalyze:
            print(f"Extracting and aligning dFF")
            self.meta = load_align.store_experiment_metadata(self.exp_info,load_info=self.params_file,
                                                            save_info=self.savepath)
            ws_data, ws_info = load_align.load_all_wavesurfer(
                self.lever_dir, min_gap, min_bout_duration_s=self.exp_info["min_img_duration_s"]
            )
            
            self.session_data,self.meta = load_align.extract_and_align_activity(
                self.img_dir, ws_info, ws_data, self.meta
            )
            
        else:
            print(f"📂 Loading existing data")
            if not (os.path.exists(meta_path) and os.path.exists(data_path)):
                raise FileNotFoundError(f"Missing session files in {self.savepath}")
            with open(meta_path, "r") as f:
                self.meta = json.load(f)
            with open(data_path, "rb") as f:
                self.session_data = pickle.load(f)
        print("✅ Session loaded and ready.")

