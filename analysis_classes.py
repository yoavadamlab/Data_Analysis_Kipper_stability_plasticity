################ Classes ################
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.pardir)))
from analysis_funcs_and_consts import *

class DB_manager:
    def __init__(self):
        self.schema = self._get_schema()
        self.db = self._get_DB()
        self.fov_analysis_options = self._get_single_fov_analysis_options()
        self.entire_db_analysis_options = self._get_entire_db_analysis_options()
        self.fig_generator = FigGenerator(self.schema, self.db)

    def _get_schema(self):
        return [
            EXPERIMENT_DATE, CAGE, MOUSE_NAME, SEQ,  GOOD_CELLS, CELL_TYPE, STRAIN, 
            FOV, FRAME_RATE, REMAPPING, REMOVED_LAPS, COMMENTS, VIDEO_PATH
        ]

    def _get_DB(self):
        df = pd.read_csv(paths.DB_PATH, names=self.schema, index_col=False)
        df = self._append_DB_columns(df)
        df = df.sort_values(EXPERIMENT_DATE, ascending=False)
        return df

    def _append_DB_columns(self, df):
        df = self._add_longitudinal_index(df)
        return df

    def _add_longitudinal_index(self, df):
        df[SESSIONS_COUNTER] = df.groupby([CAGE, MOUSE_NAME, FOV])[EXPERIMENT_DATE].transform('nunique')
        return df
    
    def _get_single_fov_analysis_options(self):
        """
        return list of analysis options relvenat to a single FOV
        """
        return [
            CELLS_ACTIVITY, FR_AND_SUB, LAP_FR,  
            ACTIVITY_PER_LAP, LONGITUDIAL_ANALYSIS
        ]

    def _get_entire_db_analysis_options(self):
        """
        return list of analysis options relvenat to the entire DB - 
        meaning some summary figure over all the dataset
        """
        return [FR_POPULATION]

    def display_single_FOV_figure(self, fig_name, data_record):
        experiment = Experiment(self.schema, data_record)
        fig = self.fig_generator.create_fig(fig_name, experiment)
        return fig
    
    def display_FOV_images(self, data_record, fig_name=None):
        experiment = Experiment(self.schema, data_record)
        fig = self.fig_generator.create_images_fig(fig_name, experiment)
        return fig

    def get_single_record(self,cage,mouse_name,pipeline_seq):
        single_record = self.db[(self.db[CAGE]==cage) & (self.db[MOUSE_NAME]==mouse_name) & (self.db[SEQ]==pipeline_seq)]
        return single_record

    def get_single_experiment(self,cage,mouse_name,pipeline_seq):
        single_record = self.get_single_record(cage,mouse_name,pipeline_seq)
        db_record = {}
        for i,j in zip(self.schema, single_record.values[0]):
            db_record[i] = j
        exp = Experiment(self.schema, db_record)
        return exp




class Cell:
    def __init__(self, exp, cell_num):
        self.exp = exp
        self.metadata = exp.metadata
        self.cell_num = cell_num
        self.trace = self._get_trace(exp, cell_num)
        self.spikes = self._get_spikes(exp, cell_num)
        self.name = f'{self.metadata[CAGE]}_{self.metadata[MOUSE_NAME]}_{self.metadata[SEQ]}_{cell_num}'
        self.cell_type = self.metadata[CELL_TYPE]
        self.strain = self.metadata[STRAIN]
        self.remapping = self.metadata[REMAPPING]
        self.identity = f'{self.metadata[CAGE]}_{self.metadata[MOUSE_NAME]}_{self.metadata[FOV]}_{cell_num}'
        self.intseq = int(self.metadata[SEQ].split('_')[0])
        self.session_num = self._calc_session_num()
        self.times_imaged = self._times_imaged()
        self.trial_num = self._calc_trial_num()
        
    def _get_trace(self, exp, cell_num):
        traces = exp.get_traces()
        cell_idx = exp.get_cell_idx(cell_num)
        return traces[cell_idx]
       
    def _get_spikes(self, exp, cell_num):
        spikes_timming = exp.get_spikes_timming()
        cell_idx = exp.get_cell_idx(cell_num)
        return spikes_timming[cell_idx]
    
    def _get_obj_name_suffix(self, per_exp_data):
        if per_exp_data:
            name = self.exp.name
        else:
            name = self.name
        return name
    
    def remapping(self):
        return self.remapping
    
    def _calc_session_num(cell): #calculate how many times this specific cell was previously imaged, separately for remapping and non-remapping experiments
        db_conditions = {}
        db_temp = DB_Analysis().get_analysis_records(db_conditions).copy()
        db_temp[SEQ]=db_temp.apply(lambda row: int(row[SEQ].split('_')[0]),axis=1) #convert seq to int for sortin, maybe should be done generally in DB manager, instead of locally here?
        try:
            cage_name,mouse_name,fov,cell_num=cell.identity.split('_')
            if cell.remapping==False:
                prev_exps_df=db_temp[(db_temp[CAGE]==cage_name) & (db_temp[MOUSE_NAME]==mouse_name) & (db_temp[FOV]==fov) & (db_temp[SEQ]<cell.intseq) & (db_temp[REMAPPING]==False)]
            elif cell.remapping==True:
                prev_exps_df=db_temp[(db_temp[CAGE]==cage_name) & (db_temp[MOUSE_NAME]==mouse_name) & (db_temp[FOV]==fov) & (db_temp[SEQ]<cell.intseq) & (db_temp[REMAPPING]==True)]
            #count the number of times cell.cell_num appears in the column of good cells, in all the previous experiments (depends on remapping)
            session_num=prev_exps_df[GOOD_CELLS].str.count(cell_num).sum()
        except:
            session_num=0   # for cases where the cell name is not in the correct format, need fixing
        return session_num

    def _times_imaged(cell): #calculate how many times this specific cell was previously imaged, does not take into account remapping
        db_conditions = {}
        db_temp = DB_Analysis().get_analysis_records(db_conditions).copy()
        db_temp[SEQ]=db_temp.apply(lambda row: int(row[SEQ].split('_')[0]),axis=1) #convert seq to int for sortin, maybe should be done generally in DB manager, instead of locally here?
        try:
            cage_name,mouse_name,fov,cell_num=cell.identity.split('_')
            prev_exps_df=db_temp[(db_temp[CAGE]==cage_name) & (db_temp[MOUSE_NAME]==mouse_name) & (db_temp[FOV]==fov) & (db_temp[SEQ]<cell.intseq)]
            #count the number of times cell.cell_num appears in the column of good cells, in all the previous experiments (without taking remapping into acount)
            times_imaged=prev_exps_df[GOOD_CELLS].str.count(cell_num).sum()
        except:
            times_imaged=0   # for cases where the cell name is not in the correct format, need fixing
        return times_imaged


    def _calc_trial_num(cell): #calculate how many times this mouse went throguh imaging, separately for remapping and non-remapping experiments
        db_conditions = {}
        db_temp = DB_Analysis().get_analysis_records(db_conditions).copy()
        db_temp[SEQ]=db_temp.apply(lambda row: int(row[SEQ].split('_')[0]),axis=1) #convert seq to int for sortin, maybe should be done generally in DB manager, instead of locally here?
        cage_name,mouse_name=cell.identity.split('_')[0],cell.identity.split('_')[1]
        
        if cell.remapping==False:
            prev_exps_df=db_temp[(db_temp[CAGE]==cage_name) & (db_temp[MOUSE_NAME]==mouse_name) & (db_temp[SEQ]<cell.intseq) & (db_temp[REMAPPING]==False)]
        elif cell.remapping==True:
            prev_exps_df=db_temp[(db_temp[CAGE]==cage_name) & (db_temp[MOUSE_NAME]==mouse_name) & (db_temp[SEQ]<cell.intseq) & (db_temp[REMAPPING]==True)]
        return len(prev_exps_df)

    def calc_trial_num_within_same_day(cell): #calculate how many times this mouse went throguh imaging, that day, separately for remapping and non-remapping experiments
        db_conditions = {}
        db_temp = DB_Analysis().get_analysis_records(db_conditions).copy()
        db_temp[SEQ]=db_temp.apply(lambda row: int(row[SEQ].split('_')[0]),axis=1) #convert seq to int for sortin, maybe should be done generally in DB manager, instead of locally here?
        cage_name,mouse_name=cell.identity.split('_')[0],cell.identity.split('_')[1]
        
        if cell.remapping==False:
            prev_exps_df=db_temp[(db_temp[CAGE]==cage_name) & (db_temp[MOUSE_NAME]==mouse_name) & (db_temp[SEQ]<cell.intseq) & (db_temp[REMAPPING]==False) & (db_temp[EXPERIMENT_DATE]==cell.metadata[EXPERIMENT_DATE])]
        elif cell.remapping==True:
            prev_exps_df=db_temp[(db_temp[CAGE]==cage_name) & (db_temp[MOUSE_NAME]==mouse_name) & (db_temp[SEQ]<cell.intseq) & (db_temp[REMAPPING]==True) & (db_temp[EXPERIMENT_DATE]==cell.metadata[EXPERIMENT_DATE])]
        return len(prev_exps_df)
    

    def get_subthreshold_trace(self,df=None):
        if df is None:
            df=self.exp.data
        subthreshold_traces = self.exp.get_traces_without_spikes(df)
        cell_idx = self.exp.get_cell_idx(self.cell_num)
        return subthreshold_traces[cell_idx]
    
    
    def get_power_spectrum(self,lower_bound=0,upper_bound=50,smoothing_factor=50,remapping=False,partial_frames=[],df=None): # partial_frames is a list of frames, e.g [0,1000]
        if df is None:
            df=self.exp.data
        subthreshold_trace=self.get_subthreshold_trace(df)
        if (len(partial_frames)>0) & (remapping==False):
            subthreshold_trace=subthreshold_trace[partial_frames[0]:partial_frames[1]] 
        psds, freqs = psd(subthreshold_trace, self.metadata[FRAME_RATE],upper_bound) # y axis is power, x axis is frequency
        psds = psds[np.where(freqs>lower_bound)]# remove low frequencies
        freqs = freqs[np.where(freqs>lower_bound)]# remove low frequencies
        psds=smoothdata(np.real(psds),smoothing_factor)

        if remapping:
            if self.metadata[REMAPPING]:
                switch_frame = self.exp.data[self.exp.data['current_World'] == 3].index[0]
            else:
                switch_frame = int(len(subthreshold_trace)/2)
            subthreshold_trace_fam=subthreshold_trace[:switch_frame]
            subthreshold_trace_nov=subthreshold_trace[switch_frame:]
            if len(partial_frames)>0:
                subthreshold_trace_fam=subthreshold_trace_fam[partial_frames[0]:partial_frames[1]] 
                subthreshold_trace_nov=subthreshold_trace_nov[partial_frames[2]:partial_frames[3]] # in remappnig, partial_frames is a list with len 4, [start_fam, end_fam, start_nov, end_nov]
            psds_fam, freqs_fam = psd(subthreshold_trace_fam, self.metadata[FRAME_RATE],upper_bound)
            psds_nov, freqs_nov = psd(subthreshold_trace_nov, self.metadata[FRAME_RATE],upper_bound)
            psds_fam = psds_fam[np.where(freqs_fam>lower_bound)]
            freqs_fam = freqs_fam[np.where(freqs_fam>lower_bound)]
            psds_nov = psds_nov[np.where(freqs_nov>lower_bound)]
            freqs_nov = freqs_nov[np.where(freqs_nov>lower_bound)]
            psds_fam=smoothdata(np.real(psds_fam),smoothing_factor)
            psds_nov=smoothdata(np.real(psds_nov),smoothing_factor)
            return psds, freqs, psds_fam, freqs_fam, psds_nov, freqs_nov

        return psds, freqs
    
    def get_power_spectrum_normed_by_sh(self,lower_bound=0,upper_bound=50,smoothing_factor=50,remapping=False,partial_frames=[],df=None): # partial_frames is a list of frames, e.g [0,1000]
        if df is None:
            df=self.exp.data
        subthreshold_trace=self.get_subthreshold_trace(df)
        df=get_avg_spike_height_per_time(cell=self, df=df) # get the average spike height per time, to normalize the power spectrum
        subthreshold_trace=subthreshold_trace/df['mean_spike_height'] # normalize the subthreshold trace by the average spike height
        subthreshold_trace=np.array(subthreshold_trace) # convert to numpy array, in case it is a pandas series
        if (len(partial_frames)>0) & (remapping==False):
            subthreshold_trace=subthreshold_trace[partial_frames[0]:partial_frames[1]] 
        psds, freqs = psd(subthreshold_trace, self.metadata[FRAME_RATE],upper_bound) # y axis is power, x axis is frequency
        psds = psds[np.where(freqs>lower_bound)]# remove low frequencies
        freqs = freqs[np.where(freqs>lower_bound)]# remove low frequencies
        psds=smoothdata(np.real(psds),smoothing_factor)

        if remapping:
            if self.metadata[REMAPPING]:
                switch_frame = self.exp.data[self.exp.data['current_World'] == 3].index[0]
            else:
                switch_frame = int(len(subthreshold_trace)/2)
            subthreshold_trace_fam=subthreshold_trace[:switch_frame]
            subthreshold_trace_nov=subthreshold_trace[switch_frame:]
            if len(partial_frames)>0:
                subthreshold_trace_fam=subthreshold_trace_fam[partial_frames[0]:partial_frames[1]] 
                subthreshold_trace_nov=subthreshold_trace_nov[partial_frames[2]:partial_frames[3]] # in remappnig, partial_frames is a list with len 4, [start_fam, end_fam, start_nov, end_nov]
            psds_fam, freqs_fam = psd(subthreshold_trace_fam, self.metadata[FRAME_RATE],upper_bound)
            psds_nov, freqs_nov = psd(subthreshold_trace_nov, self.metadata[FRAME_RATE],upper_bound)
            psds_fam = psds_fam[np.where(freqs_fam>lower_bound)]
            freqs_fam = freqs_fam[np.where(freqs_fam>lower_bound)]
            psds_nov = psds_nov[np.where(freqs_nov>lower_bound)]
            freqs_nov = freqs_nov[np.where(freqs_nov>lower_bound)]
            psds_fam=smoothdata(np.real(psds_fam),smoothing_factor)
            psds_nov=smoothdata(np.real(psds_nov),smoothing_factor)
            return psds, freqs, psds_fam, freqs_fam, psds_nov, freqs_nov

        return psds, freqs
    
    # def get_fam_and_novel_df(self,df=None): # if no remapping, returns the first and second half of the experiment, if remapping, returns the familiar and novel parts
    #     if df is None:
    #         df=self.exp.data                               
    #     if self.metadata[REMAPPING]:
    #         fam_df = df[df['current_World'] == 1]
    #         nov_df = df[df['current_World'] == 3]
    #         return fam_df, nov_df
    #     else:
    #         df=df[df['current_World'] == 1]                                   # This function also removes black screen frames! but we decided to only remove them once 
    #         max_lap=np.max(df['lap_counter'])                                 #we bin the data on space (within 'bin_data_by_position' function) 
    #         #divide the df into two hlaves, first half and second half, based on lap number
    #         first_half=df[df['lap_counter']<=max_lap/2]
    #         second_half=df[df['lap_counter']>max_lap/2]
    #         return first_half, second_half
    
    def get_fam_and_novel_df(self,df=None): # if no remapping, returns the first and second half of the experiment, if remapping, returns the familiar and novel parts
        if df is None:
            df = self.exp.data
        if self.metadata[REMAPPING]:
            novel_start = df[df[consts.WORLD] == 3].index[0]       
            df_fam = df.iloc[:novel_start]                          
            df_fam = df_fam.reset_index(drop=True)                              #27.6 new get_fam_and_novel_df function - without black screen removal. prbobaly will
            df_nov = df.iloc[novel_start:]                                      # effect several results, but we will see.
            df_nov = df_nov.reset_index(drop=True)
            return df_fam, df_nov
        else:
            max_lap = np.max(df[consts.LAP_COUNTER])
            first_half = df[df[consts.LAP_COUNTER] <= max_lap/2]
            second_half = df[df[consts.LAP_COUNTER] > max_lap/2]
            return first_half, second_half
        
    def get_fam_and_novel_df_partial(self,laps_to_consider,df=None): 
        #take the last n laps of the familiar environment and the first n laps of the novel environment
        if df is None:
            df = self.exp.data
        if self.metadata[REMAPPING]:
            novel_start = df[df[consts.WORLD] == 3].index[0]       
            df_fam = df.iloc[:novel_start]                          
            df_fam = df_fam.reset_index(drop=True)                              
            df_nov = df.iloc[novel_start:]                                      
            df_nov = df_nov.reset_index(drop=True)
        else:
            max_lap = np.max(df[consts.LAP_COUNTER])
            df_fam = df[df[consts.LAP_COUNTER] <= max_lap/2]
            df_fam = df_fam.reset_index(drop=True)
            df_nov = df[df[consts.LAP_COUNTER] > max_lap/2]
            df_nov = df_nov.reset_index(drop=True)
        
        df_fam=df_fam[df_fam['lap_counter']>np.max(df_fam['lap_counter'])-laps_to_consider]
        df_nov['lap_counter']=df_nov['lap_counter']-np.min(df_nov['lap_counter'])+1
        df_nov=df_nov[df_nov['lap_counter']<=laps_to_consider]
        return df_fam, df_nov
    
    def get_fam_and_novel_df_no_reset_index(self,df=None): 
        if df is None:
            df = self.exp.data
        if self.metadata[REMAPPING]:
            novel_start = df[df[consts.WORLD] == 3].index[0]       
            df_fam = df.iloc[:novel_start]                          
            df_nov = df.iloc[novel_start:]                                      
        else:
            max_lap = np.max(df[consts.LAP_COUNTER])
            df_fam = df[df[consts.LAP_COUNTER] <= max_lap/2]
            df_nov = df[df[consts.LAP_COUNTER] > max_lap/2]
        
        return df_fam, df_nov
    
    def get_fam_and_novel_df_no_reset_index_partial(self,laps_to_consider,df=None):
        if df is None:
            df = self.exp.data
        if self.metadata[REMAPPING]:
            novel_start = df[df[consts.WORLD] == 3].index[0]       
            df_fam = df.iloc[:novel_start]                          
            df_nov = df.iloc[novel_start:]                                      
        else:
            max_lap = np.max(df[consts.LAP_COUNTER])
            df_fam = df[df[consts.LAP_COUNTER] <= max_lap/2]
            df_nov = df[df[consts.LAP_COUNTER] > max_lap/2]
        
        df_fam=df_fam[df_fam['lap_counter']>np.max(df_fam['lap_counter'])-laps_to_consider]
        df_nov['lap_counter']=df_nov['lap_counter']-np.min(df_nov['lap_counter'])+1
        df_nov=df_nov[df_nov['lap_counter']<=laps_to_consider]
        return df_fam, df_nov

        

    
        
############### Numpy data ############################ 
# For handling matrices that are saved as numpy array, either per experiment or per cell (e.g speed_mat vs. FR_mat)

    def is_npy_data_exists(self, dirname, per_exp_data=False): #the bool flag tells us is this data is per experiment or per cell (e.g speed_mat vs. FR_mat)
        files = glob.glob(dirname + '\*.npy')
        exists = False
        # check if a file starting with self.name exists in the folder
        name = self._get_obj_name_suffix(per_exp_data)
        for file in files:
            filename = os.path.basename(file)
            if filename.startswith(name):
                exists = True                
                break
        return exists
    #check if the data exists in the folder, withot a for loop. didn't help much in reducing the time! why??
    def is_npy_data_exists2(self, dirname, per_exp_data=False):
        name = self._get_obj_name_suffix(per_exp_data)
        if self.metadata[REMAPPING]:
            fam_exists = os.path.exists(os.path.join(dirname, name + '_fam.npy'))
            nov_exists = os.path.exists(os.path.join(dirname, name + '_nov.npy'))
            return fam_exists and nov_exists
        else:
            return os.path.exists(os.path.join(dirname, name + '.npy'))
        
    def get_npy_object(self, dirname, per_exp_data=False):
        name = self._get_obj_name_suffix(per_exp_data)
        if self.metadata[REMAPPING]:
            fam = np.load(os.path.join(dirname, name + '_fam.npy'), allow_pickle=True)
            nov = np.load(os.path.join(dirname, name + '_nov.npy'), allow_pickle=True)
            return [fam, nov]
        else:
            fr_mat = np.load(os.path.join(dirname, name + '.npy'), allow_pickle=True)
            return [fr_mat]
                       
    def save_npy_object(self, obj_lst, dirname, per_exp_data=False):
        name = self._get_obj_name_suffix(per_exp_data)
        if self.metadata[REMAPPING]:
            np.save(os.path.join(dirname, name + '_fam.npy'), obj_lst[0])
            np.save(os.path.join(dirname, name + '_nov.npy'), obj_lst[1])
        else:
            np.save(os.path.join(dirname, name + '.npy'), obj_lst[0])  

################ End of Numpy data ####################


################# Json data #########################
# For handling analysis results that are saved as json file (e.g. spatial_info)

    def is_test_result_exists(self, test_name): #the bool flag tells us is this data is per experiment or per cell (e.g speed_mat vs. FR_mat)
        json_path = os.path.join(paths.ANALYSIS_RESULTS_DIR, test_name + '.json')
        if os.path.exists(json_path):
            with open(json_path) as json_file:
                results_dict = json.load(json_file)
            if self.name in results_dict.keys():
                return True 
        else:
            return False
        
        
    def get_data_from_json(self, test_name):
        json_path = os.path.join(paths.ANALYSIS_RESULTS_DIR, test_name + '.json')
        with open(json_path) as json_file:
            results_dict = json.load(json_file)
        return results_dict[self.name]   
    
                    
################ End of Json data ####################
    

    ################# Pkl data #########################
    def does_cell_result_exist_pkl(self, test_name):
        if DB_Analysis().does_test_result_exist_pkl(test_name):
            result_dict=DB_Analysis().get_data_from_pkl(test_name)
            if self.name in result_dict.keys():
                return True
            else:
                return False
        else:
            print('pkl file does not exist for entire test - ' + str(test_name))
            return False
    
    def get_cell_data_from_pkl(self, test_name):
        result_dict=DB_Analysis().get_data_from_pkl(test_name)
        return result_dict[self.name]
    
    
################# all np arrays creation #########################

    def calc_FR_matrix(self): 
        FR_matrices = self.exp.get_fr_matrices()
        if len(FR_matrices) == 1:
            cell_FR_matrix = FR_matrices[0][self.cell_num]
            # cell_FR_matrix =cell_FR_matrix[:,:-1] # remove last bin
            return [cell_FR_matrix]
        else:
            fam_FR_matrix = FR_matrices[0][self.cell_num]
            nov_FR_matrix = FR_matrices[1][self.cell_num]
            # fam_FR_matrix = fam_FR_matrix[:,:-1] # remove last bin
            # nov_FR_matrix = nov_FR_matrix[:,:-1] # remove last bin
            return [fam_FR_matrix, nov_FR_matrix]

    def get_FR_matrix(self):
        if self.is_npy_data_exists2(paths.FR_MATS_DIR):
            return self.get_npy_object(paths.FR_MATS_DIR)
        else: # calc and save
            print('FR_mat does not exist for ' + self.name + '. creating...')
            FR_mats = self.calc_FR_matrix() # returns a list of 1 or 2 matrices (depends on remapping)
            self.save_npy_object(FR_mats, paths.FR_MATS_DIR)
            return FR_mats
        
    def get_FR_matrix_new(self):
        if self.does_cell_result_exist_pkl(Cell.calc_FR_matrix):
            return self.get_cell_data_from_pkl(Cell.calc_FR_matrix)
        else:
            print('FR_mat does not exist for ' + self.name + '. creating...')
            FR_mats = self.calc_FR_matrix()
            return FR_mats
    
    def calc_time_matrix(self):
        time_mats = self.exp.calc_time_matrix_of_exp()
        return time_mats
        
    def get_time_matrix(self):
        if self.is_npy_data_exists(paths.TIME_MATS_DIR, per_exp_data=True):
            return self.get_npy_object(paths.TIME_MATS_DIR, per_exp_data=True)
        else: # calc and save
            print('time_mat does not exist for ' + self.name + '. creating...')
            time_mats = self.calc_time_matrix() # returns a list of 1 or 2 matrices (depends on remapping)
            self.save_npy_object(time_mats, paths.TIME_MATS_DIR, per_exp_data=True)
            return time_mats
        
    def get_time_matrix_new(self):
        if self.does_cell_result_exist_pkl(Cell.calc_time_matrix):
            return self.get_cell_data_from_pkl(Cell.calc_time_matrix)
        else:
            print('time_mat does not exist for ' + self.name + '. creating...')
            time_mats = self.calc_time_matrix()
            return time_mats

    def calc_speed_matrix(self):
        speed_mats= self.exp.calc_speed_matrix_of_exp()
        return speed_mats
        
    def get_speed_matrix(self):
        if self.is_npy_data_exists(paths.SPEED_MATS_DIR, per_exp_data=True):
            return self.get_npy_object(paths.SPEED_MATS_DIR, per_exp_data=True)
        else:
            print('speed_mat does not exist for ' + self.name + '. creating...')
            speed_mats = self.calc_speed_matrix()
            self.save_npy_object(speed_mats, paths.SPEED_MATS_DIR, per_exp_data=True)
            return speed_mats
        
    def get_speed_matrix_new(self):
        if self.does_cell_result_exist_pkl(Cell.calc_speed_matrix):
            return self.get_cell_data_from_pkl(Cell.calc_speed_matrix)
        else:
            print('speed_mat does not exist for ' + self.name + '. creating...')
            speed_mats = self.calc_speed_matrix()
            return speed_mats
        
################ End of all np arrays creation ####################

################# General Stats #########################

    def calc_mean_FR(self):
        FR_mats = self.get_FR_matrix()
        if len(FR_mats) == 1:
            mean_FR = np.nanmean(FR_mats[0])*1000 # convert to Hz
            return [mean_FR]
        else:
            fam_mean_FR = np.nanmean(FR_mats[0])*1000 # convert to Hz
            nov_mean_FR = np.nanmean(FR_mats[1])*1000 # convert to Hz
            return [fam_mean_FR, nov_mean_FR]
        
    def calc_mean_speed(self):
        speed_mats = self.get_speed_matrix()
        if len(speed_mats) == 1:
            mean_speed = np.nanmean(speed_mats[0])
            return [mean_speed]
        else:
            fam_mean_speed = np.nanmean(speed_mats[0])
            nov_mean_speed = np.nanmean(speed_mats[1])
            return [fam_mean_speed, nov_mean_speed]
        
    def get_spike_heights(self):
        spike_times=self.spikes
        spike_heights=self.trace[spike_times]
        return spike_heights
    
        
    def calc_SNR(self): #need to copy from Gal's analysis
        spike_heights=self.get_spike_heights()
        signal=np.mean(spike_heights) #for this calculation signal is the average spike height 
        noise=np.std(self.get_subthreshold_trace())
        SNR=signal/noise
        return SNR


    def calc_local_SNR(self):
        spikes=self.spikes
        trace=self.trace
        sub_thresh = self.get_subthreshold_trace()
        std = np.std(sub_thresh)
        all_SBRS = []
        for s in spikes:
            try:
                mini =  min(trace[s-3:s])
                sbr = ( trace[s] - mini)/std
                all_SBRS.append(sbr)
            except:
                continue
        m_SBR = np.mean(all_SBRS)
        return m_SBR
      
    def calc_isi(self):
        if self.metadata[REMAPPING]:
            remapping_frame = self.exp.data[self.exp.data['current_World'] == 3].index[0]
            diff_in_frames_fam=np.diff(self.spikes[self.spikes < remapping_frame])
            diff_in_frames_nov=np.diff(self.spikes[self.spikes > remapping_frame])
            isi_in_ms_fam=(diff_in_frames_fam/self.metadata['frame_rate'])*1000
            isi_in_ms_nov=(diff_in_frames_nov/self.metadata['frame_rate'])*1000
            return [isi_in_ms_fam, isi_in_ms_nov]
        else:
            diff_in_frames=np.diff(self.spikes)
            isi_in_ms=(diff_in_frames/self.metadata['frame_rate'])*1000
            return [isi_in_ms]

    def calc_isi_splitted(self):
        if self.metadata[REMAPPING]:
            remapping_frame = self.exp.data[self.exp.data['current_World'] == 3].index[0]
            diff_in_frames_fam=np.diff(self.spikes[self.spikes < remapping_frame])
            diff_in_frames_nov=np.diff(self.spikes[self.spikes > remapping_frame])
            isi_in_ms_fam=(diff_in_frames_fam/self.metadata['frame_rate'])*1000
            isi_in_ms_nov=(diff_in_frames_nov/self.metadata['frame_rate'])*1000
        else:
            middle_frame = int(len(self.exp.data)/2)
            diff_in_frames_fam=np.diff(self.spikes[:middle_frame])
            diff_in_frames_nov=np.diff(self.spikes[middle_frame:])
            isi_in_ms_fam=(diff_in_frames_fam/self.metadata['frame_rate'])*1000
            isi_in_ms_nov=(diff_in_frames_nov/self.metadata['frame_rate'])*1000
        return [isi_in_ms_fam, isi_in_ms_nov]
        
    def calc_isi_on_df(self,df=None,isi_threshold_ms=50):
        if df is None:
            df=self.exp.data
        spikes_col = df[SPIKES_PREFIX + str(self.cell_num)]
        spikes_times=np.where(spikes_col==1)[0]
        diff_in_frames=np.diff(spikes_times)
        isi_in_ms=(diff_in_frames/self.metadata['frame_rate'])*1000
        isi_in_ms=isi_in_ms[isi_in_ms<isi_threshold_ms]
        return isi_in_ms
    
    def calc_isi_on_df_frames(self,df=None,isi_threshold_frames=25):
        if df is None:
            df=self.exp.data
        spikes_col = df[SPIKES_PREFIX + str(self.cell_num)]
        spikes_times=np.where(spikes_col==1)[0]
        diff_in_frames=np.diff(spikes_times)
        isi_in_frames=diff_in_frames[diff_in_frames<isi_threshold_frames]
        return isi_in_frames

    
    def calc_median_isi(self):
        if self.metadata[REMAPPING]:
            isi_in_ms=self.calc_isi()
            median_isi_fam = np.median(isi_in_ms[0])
            median_isi_nov = np.median(isi_in_ms[1])
            return [median_isi_fam, median_isi_nov]
        else:
            median_isi = np.median(self.calc_isi()[0])
            return [median_isi]

    def calc_median_isi_splitted(self):
        median_isi_fam = np.median(self.calc_isi_splitted()[0])
        median_isi_nov = np.median(self.calc_isi_splitted()[1])
        return [median_isi_fam, median_isi_nov]

        
    def calc_cv_isi(self):
        if self.metadata[REMAPPING]:
            isi_in_ms=self.calc_isi()
            cv_isi_fam = np.std(isi_in_ms[0])/np.mean(isi_in_ms[0])
            cv_isi_nov = np.std(isi_in_ms[1])/np.mean(isi_in_ms[1])
            return [cv_isi_fam, cv_isi_nov]
        else:
            cv_isi = np.std(self.calc_isi()[0])/np.mean(self.calc_isi()[0])
            return [cv_isi]
    
    def calc_fano_factor(self):
        FR_mats = self.get_FR_matrix()
        if len(FR_mats) == 1:
            fano_factor = np.nanvar(FR_mats[0])/np.nanmean(FR_mats[0])
            return [fano_factor]
        else:
            fam_fano_factor = np.nanvar(FR_mats[0])/np.nanmean(FR_mats[0])
            nov_fano_factor = np.nanvar(FR_mats[1])/np.nanmean(FR_mats[1])
            return [fam_fano_factor, nov_fano_factor]
    
    def calc_burst_index(self):
            isi_in_ms = self.calc_isi()
            ten_ms_in_frames = 10*self.metadata['frame_rate']/1000
            if len(isi_in_ms) == 1:
                burst_index = np.sum(isi_in_ms[0] < ten_ms_in_frames)/len(isi_in_ms[0])
                return [burst_index]
            else:
                fam_burst_index = np.sum(isi_in_ms[0] < ten_ms_in_frames)/len(isi_in_ms[0])
                nov_burst_index = np.sum(isi_in_ms[1] < ten_ms_in_frames)/len(isi_in_ms[1])
                return [fam_burst_index, nov_burst_index]
    
    def calc_burst_index_novel_comparison(self):
        #compare the first and last third of the novel environment
        isi_in_ms = self.calc_isi()
        ten_ms_in_frames = 10*self.metadata['frame_rate']/1000
        if len(isi_in_ms) > 1:
            first_third = isi_in_ms[1][:int(len(isi_in_ms[1])/3)]
            last_third = isi_in_ms[1][int(len(isi_in_ms[1])*2/3):]
            burst_index_first_third = np.sum(first_third < ten_ms_in_frames)/len(first_third)
            burst_index_last_third = np.sum(last_third < ten_ms_in_frames)/len(last_third)

            return [burst_index_first_third, burst_index_last_third]
        
    def calc_burst_index_familiar_control(self):
        #compare the first and last sixth of the familiar environment, in experiments with no remapping
        isi_in_ms = self.calc_isi()
        ten_ms_in_frames = 10*self.metadata['frame_rate']/1000
        if len(isi_in_ms) == 1:
            first_sixth = isi_in_ms[0][:int(len(isi_in_ms[0])/6)]
            last_sixth = isi_in_ms[0][int(len(isi_in_ms[0])*5/6):]
            burst_index_first_sixth = np.sum(first_sixth < ten_ms_in_frames)/len(first_sixth)
            burst_index_last_sixth = np.sum(last_sixth < ten_ms_in_frames)/len(last_sixth)

            return [burst_index_first_sixth, burst_index_last_sixth]
    

    
    def calc_burst_index_familiar_comparison(self):
        #compare the first and last third of the novel environment
        isi_in_ms = self.calc_isi()
        ten_ms_in_frames = 10*self.metadata['frame_rate']/1000
        if len(isi_in_ms) > 1:
            first_third = isi_in_ms[0][:int(len(isi_in_ms[0])/3)]
            last_third = isi_in_ms[0][int(len(isi_in_ms[0])*2/3):]
            burst_index_first_third = np.sum(first_third < ten_ms_in_frames)/len(first_third)
            burst_index_last_third = np.sum(last_third < ten_ms_in_frames)/len(last_third)

            return [burst_index_first_third, burst_index_last_third]
        
    

    
    #Gal's functions for detecting Bursts in CKIIs:

    def get_spikes_amp(self):
        spks_amps = []
        for s in self.spikes:
            try:
                mini =  min(self.trace[s-3:s])
                amp = self.trace[s] - mini
                spks_amps.append(amp)
            except:
                continue
        return spks_amps
        
    def calculate_adp(self, spike, adp_window):
        mean_vm_before = np.mean(self.trace[spike - 5 - adp_window : spike - 4])
        mean_vm_after = np.mean(self.trace[spike + 5 : spike + 6 + adp_window ])
        adp = ((mean_vm_after - mean_vm_before) / np.mean(self.get_spikes_amp())) * 100
        return adp

    def put_sb(self, window=14):
        put_sb = {}
        fr_weight = 1 / (self.metadata[FRAME_RATE] / 1000)
        spikes = self.spikes
        cur_spike = spikes[0]
        for i in range(len(spikes)):
                if 1 <= spikes[i] - spikes[i-1] <= round (window * fr_weight):
                    if cur_spike not in put_sb.keys():
                        put_sb[cur_spike] = [cur_spike,spikes[i]]
                    else:
                        put_sb[cur_spike].append(spikes[i])  
                elif spikes[i] - spikes[i-1] > window * fr_weight:
                    cur_spike = spikes[i]
        return put_sb

    def calc_SBS(self, adp_window=20, adp_lim=15):
        SBS = {}
        s_in_bursts = []
        PSB = self.put_sb()
        for s in PSB.keys():
            if self.calculate_adp(s, adp_window) >= adp_lim:
                SBS[s] = PSB[s]
        for burst in SBS.keys():
            s_in_bursts.extend(SBS[burst])
        return SBS, s_in_bursts

    def get_SBS(self,adp_window=20, adp_lim=15):       #this takes alot of time, because the pkl files is huge... this result is way to big of an array to be saved as pkl (dictionary of all burst times + the indices of their spikes)
        if self.does_cell_result_exist_pkl(Cell.calc_SBS):     
            return self.get_cell_data_from_pkl(Cell.calc_SBS)
        else:
            SBS = self.calc_SBS(adp_window, adp_lim)
            return SBS

            
    ###################### Firing Rate per lap ######################

        def calc_fr_per_lap(self):
            FR_mats = self.get_FR_matrix()
            if len(FR_mats) == 1:
                FR_per_lap = np.nanmean(FR_mats[0], 1)
                return [FR_per_lap]
            else:
                fam_FR_per_lap = np.nanmean(FR_mats[0], 1)
                nov_FR_per_lap = np.nanmean(FR_mats[1], 1)
                return [fam_FR_per_lap, nov_FR_per_lap]
        
    
    

        
################# Spatial info ######################

    def calc_spatial_info(self):
        FR_mats = self.get_FR_matrix()
        time_mats = self.get_time_matrix()
        if len(FR_mats) == 1:
            si = self.calc_spatial_info_on_mat(FR_mats[0], time_mats[0])
            return [si]
        else:
            fam_si = self.calc_spatial_info_on_mat(FR_mats[0], time_mats[0])
            nov_si = self.calc_spatial_info_on_mat(FR_mats[1], time_mats[1])
            return [fam_si, nov_si]
        
    def calc_spatial_info_new(self):
        FR_mats = self.get_FR_matrix_new()
        time_mats = self.get_time_matrix_new()
        if len(FR_mats) == 1:
            si = self.calc_spatial_info_on_mat(FR_mats[0], time_mats[0])
            return [si]
        else:
            fam_si = self.calc_spatial_info_on_mat(FR_mats[0], time_mats[0])
            nov_si = self.calc_spatial_info_on_mat(FR_mats[1], time_mats[1])
            return [fam_si, nov_si]

    def calc_spatial_info_on_mat(self, FR_mat, time_mat):
        FR_vec = np.nanmean(FR_mat, 0)
        mean_time_vec = np.nanmean(time_mat, 0)
        pi_vec = mean_time_vec/mean_time_vec.sum() 
        r = np.dot(FR_vec,pi_vec)
        si = 0
        for ri,pi in zip(FR_vec, pi_vec):
            if ri !=0:
                si += pi*(ri/r)*np.log2(ri/r)
        return si
    
 
    def get_spatial_info(self): #the way this functions is working now, it takes a minute to run even if the data is already there, guess because we keep loading the json file?
        if self.is_test_result_exists('si'): 
            return self.get_data_from_json('si')
        else: # calc and save
            print('SI does not exist for ' + self.name + '. creating...')
            si_llst = self.calc_spatial_info() # returns a list of 1 or 2 matrices (depends on remapping) 
            return si_llst
        
    def get_spatial_info_new(self):
        if self.does_cell_result_exist_pkl(Cell.calc_spatial_info_new):
            return self.get_cell_data_from_pkl(Cell.calc_spatial_info_new)
        else:
            print('SI does not exist for ' + self.name + '. creating...')
            si_llst = self.calc_spatial_info_new() # returns a list of 1 or 2 matrices (depends on remapping) 
            return si_llst
    
    def calc_spatial_info_Zscore(self):
        FR_mats = self.get_FR_matrix()
        time_mats = self.get_time_matrix()
        if len(FR_mats) == 1:
            si_zscore = self.calc_spatial_info_Zscore_on_FR_mat(FR_mats[0], time_mats[0])
            return [si_zscore]
        else:
            fam_si_zscore = self.calc_spatial_info_Zscore_on_FR_mat(FR_mats[0], time_mats[0])
            nov_si_zscore = self.calc_spatial_info_Zscore_on_FR_mat(FR_mats[1], time_mats[1])
            return [fam_si_zscore, nov_si_zscore]
    
    def calc_spatial_info_Zscore_new(self):
        FR_mats = self.get_FR_matrix_new()
        time_mats = self.get_time_matrix_new()
        if len(FR_mats) == 1:
            si_zscore = self.calc_spatial_info_Zscore_on_FR_mat(FR_mats[0], time_mats[0])
            return [si_zscore]
        else:
            fam_si_zscore = self.calc_spatial_info_Zscore_on_FR_mat(FR_mats[0], time_mats[0])
            nov_si_zscore = self.calc_spatial_info_Zscore_on_FR_mat(FR_mats[1], time_mats[1])
            return [fam_si_zscore, nov_si_zscore]
        
    
    def calc_spatial_info_Zscore_on_FR_mat(self, FR_mat, time_mat):
        si=self.calc_spatial_info_on_mat(FR_mat, time_mat) # calc SI on original data. this has the disadvantage of not loadin SI, but this should run on single si, not on a list two si, so not for remapping
        shuffled_si_values=self.create_n_shuffled_SI_values(FR_mat, time_mat, num_of_permutations=1000 )  # to keep simulated FR vectors of every perm
        si_zscore=(si-np.mean(shuffled_si_values))/np.std(shuffled_si_values)
        
        return si_zscore

    def calc_spatial_info_and_Zscore_on_FR_mat(self, FR_mat, time_mat):
        si=self.calc_spatial_info_on_mat(FR_mat, time_mat) # calc SI on original data. this has the disadvantage of not loadin SI, but this should run on single si, not on a list two si, so not for remapping
        shuffled_si_values=self.create_n_shuffled_SI_values(FR_mat, time_mat, num_of_permutations=1000 )  # to keep simulated FR vectors of every perm
        si_zscore=(si-np.mean(shuffled_si_values))/np.std(shuffled_si_values)
        
        return [si,si_zscore]
    
    def calc_spatial_info_Pval_on_FR_mat(self, FR_mat, time_mat):    #instead of calculating Zscore, I will calculate the p value by counting the number of times 
        si=self.calc_spatial_info_on_mat(FR_mat, time_mat)          #the shuffled SI values are greater than the original SI value
        shuffled_si_values=self.create_n_shuffled_SI_values(FR_mat, time_mat, num_of_permutations=1000)  # to keep simulated FR vectors of every perm
       
        p_value = len([i for i in shuffled_si_values if i > si])/len(shuffled_si_values)
        return p_value
    
    def calc_spatial_info_Zscore_and_Pval(self,num_of_permutations=1000):
        FR_mats = self.get_FR_matrix()
        time_mats = self.get_time_matrix()
        if len(FR_mats) == 1:
            si_zscore, p_value, = self.calc_spatial_info_Zscore_and_Pval_on_FR_mat(FR_mats[0], time_mats[0], num_of_permutations)
            return [si_zscore], [p_value]
        else:
            fam_si_zscore,fam_si_pval = self.calc_spatial_info_Zscore_and_Pval_on_FR_mat(FR_mats[0], time_mats[0], num_of_permutations)
            nov_si_zscore,nov_si_pval = self.calc_spatial_info_Zscore_and_Pval_on_FR_mat(FR_mats[1], time_mats[1], num_of_permutations)
            return [fam_si_zscore, nov_si_zscore], [fam_si_pval, nov_si_pval]
    
    def calc_spatial_info_Zscore_and_Pval_on_FR_mat(self, FR_mat, time_mat, num_of_permutations):
        si=self.calc_spatial_info_on_mat(FR_mat, time_mat) # calc SI on original data. this has the disadvantage of not loadin SI, but this should run on single si, not on a list two si, so not for remapping
        shuffled_si_values=self.create_n_shuffled_SI_values(FR_mat, time_mat, num_of_permutations )  # to keep simulated FR vectors of every perm
        si_zscore=(si-np.mean(shuffled_si_values))/np.std(shuffled_si_values)
        p_value = len([i for i in shuffled_si_values if i > si])/len(shuffled_si_values)

        return si_zscore, p_value

    
    def create_n_shuffled_SI_values(self, original_FR_matrix,original_time_matrix, num_of_permutations):
        lap_number=original_FR_matrix.shape[0]
        #circshift each lap and calc SI
        shuffled_SI_values=[]
        for i in range(num_of_permutations):
            shuffled_FR_matrix=np.zeros(original_FR_matrix.shape)
            for lap in range(lap_number):
                shuffled_FR_matrix[lap,:]=np.roll(original_FR_matrix[lap,:], np.random.randint(original_FR_matrix.shape[1]))
            shuffled_SI_values.append(self.calc_spatial_info_on_mat(shuffled_FR_matrix, original_time_matrix))
        return shuffled_SI_values
    
    def get_spatial_info_Zscore(self):
        if self.does_cell_result_exist_pkl(Cell.calc_spatial_info_Zscore):
            return self.get_cell_data_from_pkl(Cell.calc_spatial_info_Zscore)
        else:
            si_zscore_llst = self.calc_spatial_info_Zscore()
            return si_zscore_llst
        
    def get_spatial_info_Zscore_new(self):
        if self.does_cell_result_exist_pkl(Cell.calc_spatial_info_Zscore_new):
            return self.get_cell_data_from_pkl(Cell.calc_spatial_info_Zscore_new)
        else:
            si_zscore_llst = self.calc_spatial_info_Zscore_new()
            return si_zscore_llst     
        
        
################# End of Spatial info ######################
        

################# in\out ratio ######################  not written yet

    def get_neighboring_indexes(self,FR_vec,threshold):
            max_index = np.argmax(FR_vec)
            neighbor_indexes = [max_index]
            start_indexes=[] #i saparated the end and start of PFs that cross from end to start. this is temporary for plotting reasons
            for i in range(max_index + 1, len(FR_vec)):
                if FR_vec[i] < threshold:
                    break
                neighbor_indexes.append(i)
            if len(FR_vec)-1 in neighbor_indexes:
                for i in range(0 , max_index):
                    if FR_vec[i] < threshold:
                        break
                    start_indexes.append(i)
                            
            for i in range(max_index - 1, -1, -1):
                if FR_vec[i] < threshold:
                    break
                neighbor_indexes.append(i)
            return np.sort(neighbor_indexes),np.sort(start_indexes)
        
    def in_out_ratio_per_cell(self,FR_vec):  
            FR_vec=smooth(FR_vec,3)
            threshold=FR_vec.mean() + 0.10 * (FR_vec.max() - FR_vec.mean())
            a,b=self.get_neighboring_indexes(FR_vec,threshold)
            norm_vec=(FR_vec-FR_vec.min())/(FR_vec.max()-FR_vec.min())#adding normalization
            if len(b)==0:
                in_field_norm_vec=norm_vec[a]
                in_field_indices=a
            else:
                in_field_norm_vec = np.concatenate((norm_vec[a],norm_vec[b]),axis=0)
                in_field_indices=np.concatenate((a,b),axis=0)
            out_field_norm_vec = [ elem for elem in norm_vec if elem not in in_field_norm_vec]
            in_fr=np.mean(in_field_norm_vec)
            out_fr=np.mean(out_field_norm_vec)
            ratio=in_fr/out_fr
            return ratio,in_field_indices,in_fr,out_fr
    
    def detect_pf_per_cell(self,FR_vec):
            ratio,in_field_indices,in_fr,out_fr=self.in_out_ratio_per_cell(FR_vec)
            return ratio,in_field_indices
        
    def detect_pf_per_cell_if_PC(self,FR_vec):
            ratio,in_field_indices,in_fr,out_fr=self.in_out_ratio_per_cell(FR_vec)
            if ratio>2.5:
                return ratio,in_field_indices
            
    
    def calc_in_out_ratio(self):
        def _calc_in_out_ratio_on_FR_mat(FR_mat):
            FR_vec = np.nanmean(FR_mat, 0)
            ratio,in_field_indices,in_fr,out_fr=self.in_out_ratio_per_cell(FR_vec)
            return ratio
        FR_mats = self.get_FR_matrix_new()
        ratios_lst=[]
        for FR_mat in FR_mats:
            ratios_lst.append(_calc_in_out_ratio_on_FR_mat(FR_mat))
        return ratios_lst

    def calc_PF_bins(self):
        def _calc_PF_bins_on_FR_mat(FR_mat):
            FR_vec = np.nanmean(FR_mat, 0)
            ratio,in_field_indices,in_fr,out_fr=self.in_out_ratio_per_cell(FR_vec)
            return in_field_indices
        # FR_mats = self.get_FR_matrix_new()
        FR_mats = self.calc_FR_matrix()
        indices_lst=[]
        for FR_mat in FR_mats:
            indices_lst.append(_calc_PF_bins_on_FR_mat(FR_mat))
        return indices_lst
    

    def get_in_out_ratio(self):
        if self.does_cell_result_exist_pkl:
            return self.get_cell_data_from_pkl(Cell.calc_in_out_ratio)
        else:
            ratio_llst = self.calc_in_out_ratio()
            return ratio_llst
    
    def get_PF_bins(self):
        if self.does_cell_result_exist_pkl:
            return self.get_cell_data_from_pkl(Cell.calc_PF_bins)
        else:
            indices_llst = self.calc_PF_bins()
            return indices_llst
        
        
    def detect_pf_per_exp(self):
            pf_dict={}
            mean_fr_df = self.get_mean_firing_rate(BINS_NUM)
            for cell in self.metadata[GOOD_CELLS]:
                FR_vec = mean_fr_df[MEAN_FR_PREFIX + str(cell)]
                FR_vec=np.array(FR_vec)
                pf_dict[cell]=self.detect_pf_per_cell(FR_vec)
            return pf_dict
        
    def get_pf_range_in_cm(self):
            pf_range_dict={}
            pf_dict = self.detect_pf_per_exp()
            cm_per_bin = self.calculate_bin_length(BINS_NUM, self.data)
            for key in pf_dict.keys():
                indxs=pf_dict[key][1]
                indx_in_cm = indxs*cm_per_bin
                pf_range_dict[key]=[indx_in_cm[1],indx_in_cm[-1]]
            return pf_range_dict
    
    def calc_reliability(self,FR_mats=None):
        rel=[]
        if FR_mats is None:
            FR_mats=self.get_FR_matrix()
        rel=self.calc_reliability_on_FR_mat(FR_mats)
        return rel
    
    def calc_reliability_on_FR_mat(self,FR_mats):
        rel=[]
        for FR_matrix in FR_mats:
            FR_vec=FR_matrix.mean(axis=0)
            detect_pf=self.detect_pf_per_cell(FR_vec)
            if detect_pf:
                cell_pf=detect_pf[1]
                mean_fr=np.mean(FR_matrix)
                std_fr=np.std(FR_matrix)
                significant_bins=np.where(FR_matrix>mean_fr+2*std_fr)
                bool_mask=np.isin(significant_bins[1], cell_pf)
                laps_with_sig_bins=significant_bins[0][bool_mask]
                laps_with_sig_bins=np.unique(laps_with_sig_bins)
                reliability=len(laps_with_sig_bins)/len(FR_matrix)*100
                rel.append(reliability)
            else:
                rel.append(0)
        return rel

    def is_place_cell(self):
        res=[]
        Z_score=self.get_spatial_info_Zscore_new()
        ratio=self.get_in_out_ratio()
        for i in range(len(Z_score)): #if remapping, Z_score & ratio are both lists with 2 elements, and we return a list of 2 bools.
            if Z_score[i]>1.69 and ratio[i]>2.5:
                res.append(True)
            else:
                res.append(False)
        return res 
    
    def calc_onset_lap(self):      #notice this fuction estimates significant bins by looknig at the entire trial, not within each lap
        def _calc_onset_lap_on_arr(arr_of_significant_laps):
            onset_lap=None   #is it a good idea for the default value to be NaN?
            for lap in arr_of_significant_laps:
                next_6_laps=np.arange(lap+1,lap+7)
                bool_array=np.isin(next_6_laps, arr_of_significant_laps)
                if np.sum(bool_array)>=3: #if there are 3 or more laps with sig. bins in the next 6 laps
                    onset_lap=lap
                    break
            return onset_lap
        FR_mat_lst=self.get_FR_matrix_new()
        PF_lst=self.get_PF_bins()
        onset_laps_lst=[]  #empty list to fill with [fam] or [fam,nov]
        for i in range(len(FR_mat_lst)): #if remapping, len(FR_mat_lst)=2
            FR_mat=FR_mat_lst[i]
            mean_fr=np.mean(FR_mat);std_fr=np.std(FR_mat)
            significant_bins=np.where(FR_mat>mean_fr+1.5*std_fr)   # find bins with FR>mean+1.5*std
            bool_mask=np.isin(significant_bins[1], PF_lst[i])   # find bins that are in the PF
            laps_with_sig_bins=significant_bins[0][bool_mask]  # find laps with sig. bins in PF
            laps_with_sig_bins=np.unique(laps_with_sig_bins)  # remove duplicates
            onset_lap=_calc_onset_lap_on_arr(laps_with_sig_bins)
            onset_laps_lst.append(onset_lap)
        return onset_laps_lst


    



################# End of in\out ratio ######################


    ################# Speed tuning ######################
    def calc_matrices_speed_corr(self):
        FR_mats = self.get_FR_matrix_new()
        speed_mats = self.get_speed_matrix_new()
        if len(FR_mats) == 1:
            speed_corr =  np.corrcoef(FR_mats[0].flatten(), speed_mats[0].flatten())[0, 1]
            return [speed_corr]
        else:
            fam_speed_corr = np.corrcoef(FR_mats[0].flatten(), speed_mats[0].flatten())[0, 1]
            nov_speed_corr = np.corrcoef(FR_mats[1].flatten(), speed_mats[1].flatten())[0, 1]
            return [fam_speed_corr, nov_speed_corr]
        
    def calc_speed_corrs(self):
        if self.remapping==False:
            speed_corr_ans=self.calc_speed_corrs_on_df()
            return [speed_corr_ans]
        else:
            fam_df, nov_df = self.exp.get_fam_and_novel_df()
            fam_speed_corr_ans = self.calc_speed_corrs_on_df(fam_df)
            nov_speed_corr_ans = self.calc_speed_corrs_on_df(nov_df)
            return [fam_speed_corr_ans, nov_speed_corr_ans]
        
    def get_speed_corrs(self):
        if self.does_cell_result_exist_pkl:
            return self.get_cell_data_from_pkl(Cell.calc_speed_corrs)
        else:
            speed_corr_llst = self.calc_speed_corrs()
            return speed_corr_llst
        

        
    def calc_speed_corrs_on_df(self,df=None):
        if df is None:
            df=self.exp.data
        # df=df[(~df['position'].between(107,128))] #removing RZ bins
        df=df[df['speed']>=20]
        spikes_col="spikes_cell_"+str(self.cell_num)
        #bin the data into 100ms bins
        df=df.reset_index(drop=True)
        df['time']=df['TS_time']-df['TS_time'][0]
        df['binned_time']=pd.cut(df['time'],bins=round(len(df)/self.metadata['frame_rate']),labels=False)
        df['spikes_per_bin']=df.groupby('binned_time')[spikes_col].transform('sum')
        df['speed_per_bin']=df.groupby('binned_time')['speed'].transform('mean') #count the number of spikes and avg speed in each bin
        df['binned_speed']=pd.cut(df['speed_per_bin'],bins=20,labels=False)  #bin the speed per bin into 20 bins
        spikes_per_bin=df.groupby('binned_speed')['spikes_per_bin'].sum()  #count the number of spikes in each speed bin
        frames_per_bin=df.groupby('binned_speed')['spikes_per_bin'].count() #count the number of frames in each bin
        #take into account only bins with more than 1% of the frames
        spikes_per_bin=spikes_per_bin[frames_per_bin>0.01*len(df)]
        frames_per_bin=frames_per_bin[frames_per_bin>0.01*len(df)]
        #divide sum of spikes by num of frames to get the firing rate in each bin
        FR_per_bin=(spikes_per_bin/frames_per_bin)
        # FR_per_bin=smooth(FR_per_bin,3)
        #turn ndarray to pandas series
        FR_per_bin=pd.Series(FR_per_bin,index=spikes_per_bin.index)
        inst_r=round(np.corrcoef(df['spikes_per_bin'],df['speed_per_bin'])[0,1],2)
        slope, intercept, r_value, p_value, std_err = stats.linregress(FR_per_bin.index,FR_per_bin)
        dict_of_ans={"inst_r":inst_r,"r_value":r_value,"p_value":p_value}
        
        return dict_of_ans
    

    ################# End of Speed tuning ######################

    ################# Vm Ramps ######################


    def get_avg_Vm_ramp_by_position(cell,bins_num=100,df=None):
        """
        This function calculates the average Vm ramp for each bin. bittner used 100 bins for 175 cm maze
        """
        if df is None:
            df = cell.exp.data
        df=cell.exp.bin_data_by_position(bins_num,df)
        trace_col_prefix = data_utils.get_trace_col_prefix(df)
        sub_trace = delete_spikes_bittner_2017(cell, df[trace_col_prefix+str(cell.cell_num)], sampling_rate=cell.metadata[FRAME_RATE])
        filtered_trace_bittner = low_pass_filter_bittner2017(sub_trace, cutoff_freq=3.0, window_size=0.2, sampling_rate=cell.metadata[FRAME_RATE])
        df["filtered_trace_bittner"] = filtered_trace_bittner
        avg_Vm_ramp_by_position = df.groupby('binned_position')['filtered_trace_bittner'].mean()
        return avg_Vm_ramp_by_position



        


        
class Experiment:
    def __init__(self, db_schema, db_record):
        self.metadata = self._extract_metadata_from_db_record(db_schema, db_record)
        self.raw_data = self._get_experiment_data()
        self.data = self.preprocessing(self.raw_data)
        self.behavior_flag = data_utils.check_behavioral_data(self.data)
        self.cells = self._create_cells()
        self.name = f'{self.metadata[CAGE]}_{self.metadata[MOUSE_NAME]}_{self.metadata[SEQ]}'

    def _create_cells(self):
        cells = {}
        for cell_num in self.metadata[GOOD_CELLS]:
            cell = Cell(self, cell_num)
            cells[cell_num] = cell
        return cells
    
    def _extract_metadata_from_db_record(self, db_schema, db_record):
        metadata = {}
        for field in db_schema:
            if field == GOOD_CELLS:
                good_cells = [n for n in ast.literal_eval(db_record[field])]
                good_cells.sort()
                metadata[field] = good_cells
                metadata[CELLS_NUM] = len(good_cells)
                continue
            if field == REMAPPING:
                val = db_record[field]
                if isinstance(val, bool):
                    metadata[field] = val
                else:
                    metadata[field] = val.lower() == "true"
                continue
            metadata[field] = db_record[field]
        return metadata
    
    def _get_experiment_data(self):
        """
        get the merged_data object.
        slice out the first 1000 rows.
        """
        data_file = os.path.join(paths.DATASET_DIR_WINDOWS, self.metadata[CAGE], self.metadata[MOUSE_NAME], self.metadata[SEQ] + '.parquet')
        df = pd.read_parquet(data_file)
        return df
    
    def preprocessing(self, df):
        # 1. cut first 1000 frames detrend traces and transform to delta F
        df = self.preprocess_traces(df)
        # for behavioral experiments:
        if data_utils.check_behavioral_data(df):
            # 1. delete bad laps 
            df["min_position"] = df.groupby(consts.LAP_COUNTER)[consts.POSITION].transform(min)  # new column for removing bad laps (with position under -30)
            df = df[df["min_position"] > -35] # virmen unit that should't be crossed
            # 3. slice out all frames before the animal started to walk
            df["chnged_position"] = df[consts.POSITION].astype(float).diff()
            start = df[df["chnged_position"] != 0].index[1]
            df = df.iloc[start:]
            df.reset_index(inplace=True, drop=True)
            # 4. slice out first and last laps                                                              
            df = self._slice_out_first_and_last_laps(df)
        return df
    
    def transform_to_df_over_f(self, trace):
        trace = trace + (2*np.abs(trace.min()))  # prevent negative values
        # F0 = trace[-1000:].mean()     # calculate F0 as the mean of the last 1000 frames 
        F0 = np.percentile(trace, 10)   # Rotem changed to 10th percentile as the F0 23.12.24
        trace = (trace - F0) / F0
        return trace
    
    def preprocess_traces(self, df):
        df = df.iloc[1000:, :]
        df.reset_index(inplace=True, drop=True)
        trace_col_prefix = data_utils.get_trace_col_prefix(df)
        for cell in self.metadata[GOOD_CELLS]:
            trace = df[trace_col_prefix + str(cell)]
            try: # try to detrend if doesnt work keep the non detrend trace
                trace = data_utils.detrend_func(np.arange(len(trace)), trace)
            except:
                print("detrend problem with cell number: ", cell)
            trace = self.transform_to_df_over_f(trace)
            df[trace_col_prefix + str(cell)] = trace
        return df
    
    def get_fr_matrices(self):
        def fr_per_lap_df_to_fr_matrices(fr_per_lap_df):
            fr_per_lap_df.sort_index(inplace=True)
            # create an empty dictionary to store the numpy arrays
            FR_mats = {}
            # iterate over the cell columns in the DataFrame
            for col_name in fr_per_lap_df.columns:
                if col_name.startswith('fr_cell_'):
                    # pivot the DataFrame on lap_counter and binned_position for the current cell column
                    pivot_df = fr_per_lap_df.pivot_table(index='lap_counter', columns='binned_position', values=col_name)
                    # convert the pivot table to a numpy array
                    FR_matrix = pivot_df.to_numpy()
                    FR_matrix = np.nan_to_num(FR_matrix)
                    # store the numpy array in the dictionary using the cell name as the key
                    cell_name = int(col_name.replace('fr_cell_', ''))
                    FR_mats[cell_name] = FR_matrix
            return FR_mats
        
        bins_num = BINS_NUM
        if self.metadata[REMAPPING]:
            fam_df, nov_df = self.get_fam_and_novel_df()
            fam_fr_per_lap_df = self.get_firing_rate_per_lap(bins_num, fam_df)
            nov_fr_per_lap_df = self.get_firing_rate_per_lap(bins_num, nov_df)
            return [fr_per_lap_df_to_fr_matrices(fam_fr_per_lap_df), fr_per_lap_df_to_fr_matrices(nov_fr_per_lap_df)]
        
        fr_per_lap_df = self.get_firing_rate_per_lap(bins_num)
        return [fr_per_lap_df_to_fr_matrices(fr_per_lap_df)]

    def get_spike_shapes(self,df=None): #return a dictionary with the spike shapes of each cell, the key is the cell number, doen't include complex spikes
        if df is None:
            df = self.data
        print(self.name)
        df = self.preprocess_traces(self.raw_data)
        traces = self.get_traces(df)
        spikes_timming = self.get_spikes_timming(df)
        spike_shapes=[]
        hundred_ms = round(self.metadata[FRAME_RATE]*0.1)
        for i, cell in enumerate(self.metadata[GOOD_CELLS]):
            spike_vecs = []
            trace = traces[i]
            spikes_times = spikes_timming[i]
            for j, spike in enumerate(spikes_times[:-1]):
                if spike-hundred_ms<0 or spike+hundred_ms>len(trace) or spikes_times[j+1]-spike<hundred_ms or spike-spikes_times[j-1]<hundred_ms: #don't include complex spikes
                    continue
                spike_vector = trace[spike-hundred_ms:spike+hundred_ms]
                norm_spike_vector = normalize(spike_vector)
                spike_vecs.append(norm_spike_vector)
            spike_vecs = np.array(spike_vecs)
            mean_spike=spike_vecs.mean(axis=0)
            spike_shapes.append(mean_spike)
        return spike_shapes

    
    def get_spatial_info_per_exp(self,bins_num=BINS_NUM):
        si_dict={} #list of spatial info per cell
        mean_fr_df = self.get_mean_firing_rate(bins_num)
        mean_time_vec= self.get_mean_time_per_bin(bins_num)
        for cell in self.metadata[GOOD_CELLS]:
            FR_vec = mean_fr_df[MEAN_FR_PREFIX + str(cell)]
            FR_vec=np.array(FR_vec)
            si = self.get_spatial_info_per_cell(FR_vec,mean_time_vec)
            si_dict[cell]=si
        return si_dict #not finished!!! what should i actually return? in which form and where do i store it?
    
    # def get_neighboring_indexes(self,FR_vec,threshold):
    #     max_index = np.argmax(FR_vec)
    #     neighbor_indexes = [max_index]
    #     start_indexes=[] #i saparated the end and start of PFs that cross from end to start. this is temporary for plotting reasons
    #     for i in range(max_index + 1, len(FR_vec)):
    #         if FR_vec[i] < threshold:
    #             break
    #         neighbor_indexes.append(i)
    #     if len(FR_vec)-1 in neighbor_indexes:
    #         for i in range(0 , max_index):
    #             if FR_vec[i] < threshold:
    #                 break
    #             start_indexes.append(i)
                        
    #     for i in range(max_index - 1, -1, -1):
    #         if FR_vec[i] < threshold:
    #             break
    #         neighbor_indexes.append(i)
    #     return np.sort(neighbor_indexes),np.sort(start_indexes)
    
    # def in_out_ratio_per_cell(self,FR_vec):   # for some reason this Ratio comes out lower than the one i have in my other code. 
    #     # i found that the reason is that in the other code FR_vec is smoothed by 3
    #     threshold=FR_vec.mean() + 0.10 * (FR_vec.max() - FR_vec.mean())
    #     a,b=self.get_neighboring_indexes(FR_vec,threshold)
    #     norm_vec=(FR_vec-FR_vec.min())/(FR_vec.max()-FR_vec.min())#adding normalization
    #     if len(b)==0:
    #         in_field_norm_vec=norm_vec[a]
    #         in_field_indices=a
    #     else:
    #         in_field_norm_vec = np.concatenate((norm_vec[a],norm_vec[b]),axis=0)
    #         in_field_indices=np.concatenate((a,b),axis=0)
    #     out_field_norm_vec = [ elem for elem in norm_vec if elem not in in_field_norm_vec]
    #     in_fr=np.mean(in_field_norm_vec)
    #     out_fr=np.mean(out_field_norm_vec)
    #     ratio=in_fr/out_fr
                                                                                            # ALL OF THESE WERE MOVED TO CELL CLASS, MAKE SURE NO ISSUES WITH THAT
        
    #     return ratio,in_field_indices
    
    # def detect_pf_per_cell(self,FR_vec):
    #     ratio,in_field_indices=self.in_out_ratio_per_cell(FR_vec)
    #     if ratio>3:
    #         return True,in_field_indices
    
    # def detect_pf_per_exp(self):
    #     pf_dict={}
    #     mean_fr_df = self.get_mean_firing_rate(BINS_NUM)
    #     for cell in self.metadata[GOOD_CELLS]:
    #         FR_vec = mean_fr_df[MEAN_FR_PREFIX + str(cell)]
    #         FR_vec=np.array(FR_vec)
    #         pf_dict[cell]=self.detect_pf_per_cell(FR_vec)
    #     return pf_dict
    
    # def get_pf_range_in_cm(self):
    #     pf_range_dict={}
    #     pf_dict = self.detect_pf_per_exp()
    #     cm_per_bin = self.calculate_bin_length(BINS_NUM, self.data)
    #     for key in pf_dict.keys():
    #         indxs=pf_dict[key][1]
    #         indx_in_cm = indxs*cm_per_bin
    #         pf_range_dict[key]=[indx_in_cm[1],indx_in_cm[-1]]
    #     return pf_range_dict

    

    def get_traces(self, df=None):
        """
        return np.array in the shape: #cells, #frames
        """
        if df is None:
            df = self.data
        trace_col_prefix = data_utils.get_trace_col_prefix(df)
        traces = np.empty((self.metadata[CELLS_NUM], len(df)))
        for i, cell in enumerate(self.metadata[GOOD_CELLS]):
            trace = df[trace_col_prefix + str(cell)]
            traces[i] = trace
        return traces

    def get_normalized_traces(self,df=None):
        """
        return np.array in the shape: #cells, #frames, with the traces normalized by the avg spike height
        """
        if df is None:
            df = self.data
        traces = self.get_traces(df)
        #create empty nd array in the size of the traces
        normalized_traces=np.zeros(traces.shape) 
        for i,cell in enumerate (self.metadata[GOOD_CELLS]):
            norm_trace=[]
            spikes=df[SPIKES_PREFIX+str(cell)]
            #avg the spike heights every 10K frames of the trace, normalize the trace by the avg spike height, and store in norm_trace
            for j in range(0,len(traces[i]),10000):
                partial_spikes=spikes[j:j+10000] 
                spikes_heights=traces[i][partial_spikes[partial_spikes==1].index]
                avg_spike_height=np.mean(spikes_heights)
                partial_norm_trace=traces[i][j:j+10000]/avg_spike_height
                norm_trace.extend(partial_norm_trace)
            normalized_traces[i]=norm_trace
        return normalized_traces

    def switch_trace_cols_to_normalized(self, df=None):
        """
        return df with the normalized traces (by spike height) instead of the original traces, but using the same column names,
        so all the function that us trace_col_prefix will work, when ran on the norm_df as input
        """
        if df is None:
            df = self.data
        norm_df=df.copy()
        norm_traces = self.get_normalized_traces(df)
        trace_col_prefix = data_utils.get_trace_col_prefix(df)
        for i, cell in enumerate(self.metadata[GOOD_CELLS]):
            norm_df[trace_col_prefix + str(cell)] = norm_traces[i]
        return norm_df
             

    def _slice_out_first_and_last_laps(self, df):
        """
        The first lap can be trimmed from slicing out the first 1000 
        frames in the preprocessing step.
        The lap can halt in the middle due to imaging time.
        Since those laps are not complete they mess the data and
        rge statistics and better to be removed.
        """
        df = df[df[consts.LAP_COUNTER] > df[consts.LAP_COUNTER].min()]
        df = df[df[consts.LAP_COUNTER] < df[consts.LAP_COUNTER].max()]
        df[consts.LAP_COUNTER] = df[consts.LAP_COUNTER] - df[consts.LAP_COUNTER].min() + 1
        df = df.reset_index(drop=True)
        return df

    def get_fam_and_novel_df(self, df=None):
        if df is None:
            df = self.data
        novel_start = df[df[consts.WORLD] == 3].index[0]       
        df_fam = df.iloc[:novel_start]                          
        df_fam = df_fam.reset_index(drop=True)
        df_nov = df.iloc[novel_start:]
        df_nov = df_nov.reset_index(drop=True)
        return df_fam, df_nov
    
    def get_fam_and_novel_df_on_exp(self,df=None): # if no remapping, returns the first and second half of the experiment, if remapping, returns the familiar and novel parts
        if df is None:
            df = self.data
        if self.metadata[REMAPPING]:
            novel_start = df[df[consts.WORLD] == 3].index[0]       
            df_fam = df.iloc[:novel_start]                          
            df_fam = df_fam.reset_index(drop=True)                              #27.6 new get_fam_and_novel_df function - without black screen removal. prbobaly will
            df_nov = df.iloc[novel_start:]                                      # effect several results, but we will see.
            df_nov = df_nov.reset_index(drop=True)
            return df_fam, df_nov
        else:
            max_lap = np.max(df[consts.LAP_COUNTER])
            first_half = df[df[consts.LAP_COUNTER] <= max_lap/2]
            second_half = df[df[consts.LAP_COUNTER] > max_lap/2]
            return first_half, second_half
    


    def bin_data_by_position(self, bins_num, df):
        df=df[df[consts.WORLD]!=5] #remove the black screen. added 29.5 by Rotem. make sure to check how this removal influence the FR_mats.
        df[BIN], _ = pd.cut(df[consts.POSITION], bins=bins_num, labels=np.arange(bins_num), include_lowest=True, retbins=True) #notice that before this removal, bin size was not stable across laps and experiments
        df[BIN] = pd.to_numeric(df[BIN])
        df.reset_index(drop=True, inplace=True) # Rotem added this line 9.7.24 to make sure the index is correct, cause it ruined the get_fam_novel_df function for cells
        return df  

    def calculate_bin_length(self, df=None, bins_num=BINS_NUM):
        if df is None:
            df = self.data
        df=self.bin_data_by_position(bins_num, df)
        mean_lap_length = int(df.groupby(consts.LAP_COUNTER)[consts.LAP_LEN_CUMSUM].max().drop_duplicates().mean() / 10)
        cm_per_bin = round(mean_lap_length / bins_num, 2)
        return cm_per_bin

    def get_firing_rate_per_lap(self, bins_num, df=None):
        if df is None:
            df = self.data
        # calculate time in bin per lap column
        df = self.bin_data_by_position(bins_num, df)
        laps_dfs = df.groupby([consts.LAP_COUNTER, BIN])
        df[TIME_IN_BIN_PER_LAP] = laps_dfs[consts.TS_TIME].transform(max) - laps_dfs[consts.TS_TIME].transform(min) + 0.0001
        # count spikes per bin per lap
        spike_col_prefix = data_utils.get_spike_col_prefix(df)
        spikes = [spike_col_prefix + str(i) for i in self.metadata[GOOD_CELLS]]
        spikes_in_bin_per_lap_cols = ["spikes_count_per_lap_per_bin_cell_" + str(i) for i in self.metadata[GOOD_CELLS]]
        df[spikes_in_bin_per_lap_cols] = laps_dfs[spikes].transform('sum')
        # keep one row only per lap per bin
        df = df[[consts.LAP_COUNTER, TIME_IN_BIN_PER_LAP, BIN,WORLD] + spikes_in_bin_per_lap_cols].drop_duplicates()
        # calculate firing rate
        fr_cols = [FR_PREFIX + str(i) for i in self.metadata[GOOD_CELLS]]
        df[fr_cols] = df[spikes_in_bin_per_lap_cols].div(df[TIME_IN_BIN_PER_LAP], axis=0)
        df = df[[consts.LAP_COUNTER, BIN,WORLD] + fr_cols].drop_duplicates()
        df = df.interpolate() # just filling the nan values with interpolation
        df.reset_index(drop=True, inplace=True)
        df.sort_values(BIN, inplace=True)
        np.where(df[fr_cols].isna(), 0, df[fr_cols] ) # rotem added this line to replace nan with 0 in raster plot, but it is not working
        return df
    
    def get_subthreshold_activity_per_lap(self, bins_num, df=None): #not finished yet
        if df is None:
            df = self.data
        df=self._slice_out_first_and_last_laps(df)
        df=self.remove_spike_activty(df)
        df = self.bin_data_by_position(bins_num, df)
        laps_dfs = df.groupby([consts.LAP_COUNTER,BIN]).mean()
        # generate subthreshold activity per bin per lap
        for cell in self.metadata[GOOD_CELLS]:
            subthresh_col = NON_SPIKED_TRACE_PREFIX + str(cell)
            subthresh_in_bin_per_lap_col = "subthresh_per_lap_per_bin_cell_" + str(cell)
            df[subthresh_in_bin_per_lap_col] = laps_dfs[subthresh_col].transform('mean')

        return 

    
    def get_mean_firing_rate(self, bins_num, df=None):
        if df is None:
            df = self.data
        df = self.get_firing_rate_per_lap(bins_num, df)
        # cols names
        fr_cols = [FR_PREFIX + str(i) for i in self.metadata[GOOD_CELLS]]
        mean_fr_cols = [MEAN_FR_PREFIX + str(i) for i in self.metadata[GOOD_CELLS]]
        sem_fr_cols = [SEM_FR_PREFIX + str(i) for i in self.metadata[GOOD_CELLS]]
        # calc fr
        df[mean_fr_cols] = df.groupby([BIN])[fr_cols].transform('mean') * 1000  # for Hz units
        laps_number = self.get_laps_number(df)
        df[sem_fr_cols] = df.groupby([BIN])[fr_cols].transform('std') * 1000 * (1/np.sqrt(laps_number))  # for Hz unit
        df = df[[BIN] +  mean_fr_cols + sem_fr_cols].drop_duplicates()
        df = df.interpolate() # just filling the nan values with interpolation
        df.reset_index(drop=True, inplace=True)
        df.sort_values(BIN, inplace=True)
        return df
    
    def get_mean_time_per_bin(self, bins_num=BINS_NUM, df=None):
        if df is None:
            df = self.data
        df = self.bin_data_by_position(bins_num, df)
        laps_dfs = df.groupby([consts.LAP_COUNTER, BIN])
        df[TIME_IN_BIN_PER_LAP] = laps_dfs[consts.TS_TIME].transform(max) - laps_dfs[consts.TS_TIME].transform(min) + 0.0001
        mean_time_per_bin=smooth(df.groupby(BIN).mean()[TIME_IN_BIN_PER_LAP],2)
        return mean_time_per_bin
    
    def remove_spike_activty(self, df,window_size=3): #window size is in frames, how many frames before and after the spike to interpolate. Rotem changed to 1 cause it looks better
        traces = self.get_traces(df)
        spikes = self.get_spikes_timming(df)
        for i, cell in enumerate(self.metadata[GOOD_CELLS]):
            trace = traces[i]
            for spike_time in spikes[i]:
                start_index = max(0, spike_time - window_size)
                end_index = min(len(trace) - 1, spike_time + window_size)
                trace[start_index:end_index+1] = np.nan
            indices = np.arange(len(trace))
            trace = np.interp(indices, indices[~np.isnan(trace)], trace[~np.isnan(trace)])
            df[NON_SPIKED_TRACE_PREFIX + str(cell)] = trace 
        return df

    def get_subthreshold_activity(self, bins_num, df=None):   #this is used for getting subthreshold activity with spatial bins for the DB plot
        if df is None:
            df = self.data
        df = self.bin_data_by_position(bins_num, df)
        df = self.remove_spike_activty(df)
        sub_activity_cols = [NON_SPIKED_TRACE_PREFIX + str(i) for i in self.metadata[GOOD_CELLS]]
        mean_sub_cols = [SUB_ACTIVITY_PREFIX + str(i) for i in self.metadata[GOOD_CELLS]]
        sem_sub_cols = [SEM_SUB_ACTIVITY_PREFIX + str(i) for i in self.metadata[GOOD_CELLS]]
        # calc fr
        df[mean_sub_cols] = df.groupby([BIN])[sub_activity_cols].transform('mean')
        laps_number = self.get_laps_number(df)
        df[sem_sub_cols] = df.groupby([BIN])[sub_activity_cols].transform('std') * (1/np.sqrt(laps_number))  
        df = df[[BIN] +  mean_sub_cols + sem_sub_cols].drop_duplicates()
        df = df.interpolate() # just filling the nan values with interpolation
        df.reset_index(drop=True, inplace=True)
        df.sort_values(BIN, inplace=True)
        return df
    

    
    def get_traces_without_spikes(self,df=None):
        if df is None:
            df = self.data
        removed_df = self.remove_spike_activty(df)
        subthresh_col_prefix = NON_SPIKED_TRACE_PREFIX
        traces = np.empty((self.metadata[CELLS_NUM], len(removed_df)))
        for i, cell in enumerate(self.metadata[GOOD_CELLS]):
            trace = removed_df[subthresh_col_prefix + str(cell)]
            traces[i] = trace
        return traces

        
    def get_spikes_timming(self, df=None):
        if df is None:
            df = self.data
        spike_col_prefix = data_utils.get_spike_col_prefix(df)
        spikes_timmimg = []
        for cell in self.metadata[GOOD_CELLS]:
            spikes_time = df[df[spike_col_prefix + str(cell)] == 1].index
            spikes_timmimg.append(spikes_time)
        return spikes_timmimg

    def get_spikes_height(self, spikes_timming, traces):
        heights = []
        for i, spikes_time in enumerate(spikes_timming):
            spikes_height = 1.01 * traces[i][spikes_time]
            heights.append(spikes_height)
        return heights
    
    def get_position(self, df=None):
        if df is None:
            df = self.data
        return df[consts.POSITION]
    
    def get_licks_timming(self, df=None):
        if df is None:
            df = self.data
        return df[df[consts.VIR_LICK] == 1].index
    
    def get_laps_number(self, df):
        laps_number = df[consts.LAP_COUNTER].nunique()
        return laps_number
    
    def get_reward_zone_bins(self, df=None, BINS=BINS_NUM):
        if df is None:
            df = self.data
        df = self.bin_data_by_position(BINS, df)
        rwd_df = df[[consts.POSITION, BIN]]
        start_rwd = 107 # hard coded fron VIRMEN file
        end_rwd = 128 
        start_rwd_bin = rwd_df[(rwd_df[consts.POSITION] >= start_rwd) & (rwd_df[consts.POSITION] < end_rwd)][BIN].min()
        end_rwd_bin = rwd_df[(rwd_df[consts.POSITION] >= start_rwd) & ( rwd_df[consts.POSITION] <= end_rwd)][BIN].max()
        return start_rwd_bin, end_rwd_bin
    
    def get_cell_idx(self, cell_num):
        for i, cell_number in enumerate(self.metadata[GOOD_CELLS]):
            if cell_number == cell_num:
                return i
            
            
    def calc_time_matrix_of_exp(self):
        def calc_time_matrix_of_df(df):
            df = self.bin_data_by_position(BINS_NUM, df)
            laps_dfs = df.groupby([consts.LAP_COUNTER, BIN])
            df[TIME_IN_BIN_PER_LAP] = laps_dfs[consts.TS_TIME].transform(max) - laps_dfs[consts.TS_TIME].transform(min) + 0.0001
            df.groupby([consts.LAP_COUNTER, BIN]).mean()[TIME_IN_BIN_PER_LAP]
            # transform into a matrix of size (laps, bins)
            time_mat=df.groupby([consts.LAP_COUNTER, BIN]).mean()[TIME_IN_BIN_PER_LAP].unstack().values
            #remove last bin from the matrix
            time_mat = time_mat[:,:-1]
            return time_mat
        
        
        if self.metadata[REMAPPING]:
            fam_df, nov_df = self.get_fam_and_novel_df()
            time_mat_fam = calc_time_matrix_of_df(fam_df)
            time_mat_nov = calc_time_matrix_of_df(nov_df)
            return [time_mat_fam, time_mat_nov]
        df = self.data
        time_mat = calc_time_matrix_of_df(df)
        return [time_mat]
    
    def calc_speed_matrix_of_exp(self):
        bins_num = BINS_NUM
        def calc_speed_matrix_of_df(df):
            df = self.bin_data_by_position(bins_num, df)
            laps_dfs = df.groupby([consts.LAP_COUNTER, BIN])
            df[SPEED_IN_BIN_PER_LAP] = laps_dfs[consts.SPEED].transform('mean')
            df.groupby([consts.LAP_COUNTER, BIN]).mean()[SPEED_IN_BIN_PER_LAP]
            # transform into a matrix of size (laps, bins)
            speed_mat=df.groupby([consts.LAP_COUNTER, BIN]).mean()[SPEED_IN_BIN_PER_LAP].unstack().values
            #remove last bin from the matrix
            speed_mat = speed_mat[:,:-1]
            return speed_mat
        
        if self.metadata[REMAPPING]:
            fam_df, nov_df = self.get_fam_and_novel_df()
            speed_mat_fam = calc_speed_matrix_of_df(fam_df)
            speed_mat_nov = calc_speed_matrix_of_df(nov_df)
            return [speed_mat_fam, speed_mat_nov]
        df = self.data
        speed_mat = calc_speed_matrix_of_df(df)
        return [speed_mat]
    
    def get_pairwise_correlation(self,params={'smoothing':5,'freqs':False},df=None):  
        if df is None:
            df = self.data
        all_corrs=[]
        df = self.remove_spike_activty(df)
        sub_traces = self.get_traces_without_spikes(df)
        if len(self.metadata[GOOD_CELLS]) == 1:
            print(f'only one cell in {self.name}')
            return np.nan                               #don't forget that this returns Nans for experiments with only one cell
        # trace_color = fig_gen.get_trace_color(self)
        # sub_trace_color = fig_gen.get_subthreshold_color()
        # spike_color = fig_gen.get_spikes_color(self)
        for i, cell in enumerate(self.metadata[GOOD_CELLS]):
            sub_trace = sub_traces[i]
            sub_trace = smooth(sub_trace, params['smoothing']) # smoothing the subthreshold in order to get rid of spikes that are not detected by the algorithm and remove noise.
            corrs=[]
            for j, cell2 in enumerate(self.metadata[GOOD_CELLS]):
                    if i != j:
                        sub_trace2 = sub_traces[j]
                        sub_trace2 = smooth(sub_trace2, params['smoothing'])
                        if params['freqs']!=False:
                            sub_trace=butter_bandpass_filter(sub_trace, params['freqs'][0], params['freqs'][1], self.metadata[FRAME_RATE], order=3)
                            sub_trace2=butter_bandpass_filter(sub_trace2, params['freqs'][0], params['freqs'][1], self.metadata[FRAME_RATE], order=3)
                        # #count nans in each trace, print the number of nans and the percentage of nans, and remove them
                        # nans1=np.count_nonzero(np.isnan(sub_trace))
                        # nans2=np.count_nonzero(np.isnan(sub_trace2))
                        # if nans1 != 0:
                        #     print(f'cell {cell} has {nans1} nans, which is {round(nans1/len(sub_trace)*100,2)}% of the trace')
                        #     print (self.name)
                        # if nans2 != 0:
                        #     print(f'cell {cell2} has {nans2} nans, which is {round(nans2/len(sub_trace2)*100,2)}% of the trace')
                        #     print (self.name)
                        # #find all the indices of nan value in either trace
                        # nan_indices1=np.argwhere(np.isnan(sub_trace))
                        # nan_indices2=np.argwhere(np.isnan(sub_trace2))
                        # indxs_to_remove=np.concatenate((nan_indices1,nan_indices2),axis=0)                               ### deleted this part for now cause it doesn't seem like we have nans in the traces
                        # indxs_to_remove=np.unique(indxs_to_remove)
                        # sub_trace=np.delete(sub_trace,indxs_to_remove)
                        # sub_trace2=np.delete(sub_trace2,indxs_to_remove)
                        corr=np.corrcoef(sub_trace,sub_trace2)[0,1]
                        corrs.append(corr)
            all_corrs.append(corrs)
        all_corrs=np.round(all_corrs,5)
        return np.unique(all_corrs) #to get rid of the duplicates
    
    def get_pairwise_correlation_with_remapping(self,params={'smoothing':5,'freqs':False},df=None,partial_laps=False): #partial laps should be list laps to include in the analysis [frst_lap_fam, last_lap_fam, frst_lap_nov, last_lap_nov]
        if df is None:
            df = self.data
        if self.metadata[REMAPPING]==True:
            fam_df, nov_df = self.get_fam_and_novel_df(df)
            if partial_laps:
                fam_df=fam_df[(fam_df[consts.LAP_COUNTER]>=partial_laps[0]) & (fam_df[consts.LAP_COUNTER]<=partial_laps[1])]
                frst_lap_nov=nov_df[consts.LAP_COUNTER].min()
                nov_df=nov_df[(nov_df[consts.LAP_COUNTER]>=frst_lap_nov+partial_laps[2]) & (nov_df[consts.LAP_COUNTER]<=frst_lap_nov+partial_laps[3])]
            fam_corrs=self.get_pairwise_correlation(params,fam_df)
            nov_corrs=self.get_pairwise_correlation(params,nov_df)
            return fam_corrs,nov_corrs
        else:
            #split the data into two halves
            fam_df1=df.iloc[:int(len(df)/2)]
            fam_df2=df.iloc[int(len(df)/2):]
            if partial_laps:
                fam_df1=fam_df1[(fam_df1[consts.LAP_COUNTER]>=partial_laps[0]) & (fam_df1[consts.LAP_COUNTER]<=partial_laps[1])]
                frst_lap_fam2=fam_df2[consts.LAP_COUNTER].min()
                fam_df2=fam_df2[(fam_df2[consts.LAP_COUNTER]>=frst_lap_fam2+partial_laps[2]) & (fam_df2[consts.LAP_COUNTER]<=frst_lap_fam2+partial_laps[3])]
            fam_corrs1=self.get_pairwise_correlation(params,fam_df1)
            fam_corrs2=self.get_pairwise_correlation(params,fam_df2)

            return fam_corrs1,fam_corrs2


    def get_pairwise_correlation_comparing_parts_in_nov(self,params={'smoothing':5,'freqs':False},df=None,laps_in_nov=3,comp_to_last_only=False): #number of first laps in nov to compare to all others
        if df is None:
            df = self.data
        if self.metadata[REMAPPING]==True:
            fam_df, nov_df = self.get_fam_and_novel_df(df)
            frst_lap_nov=nov_df[consts.LAP_COUNTER].min()
            frst_laps_nov_df=nov_df[(nov_df[consts.LAP_COUNTER]>=frst_lap_nov) & (nov_df[consts.LAP_COUNTER]<=frst_lap_nov+laps_in_nov-1)]
            last_laps_nov_df=nov_df[nov_df[consts.LAP_COUNTER]>=frst_lap_nov+laps_in_nov]
            if comp_to_last_only:
                last_laps_nov_df=last_laps_nov_df[(last_laps_nov_df[consts.LAP_COUNTER]>=last_laps_nov_df[consts.LAP_COUNTER].max()-laps_in_nov+1)]
            frst_laps_nov_corrs=self.get_pairwise_correlation(params,frst_laps_nov_df)
            last_laps_nov_corrs=self.get_pairwise_correlation(params,nov_df)
            return frst_laps_nov_corrs,last_laps_nov_corrs
        else:
            #split the data into two halves, examine only second half
            fam_df2=df.iloc[int(len(df)/2):]
            frst_lap_fam2=fam_df2[consts.LAP_COUNTER].min()
            frst_laps_fam2_df=fam_df2[(fam_df2[consts.LAP_COUNTER]>=frst_lap_fam2) & (fam_df2[consts.LAP_COUNTER]<=frst_lap_fam2+laps_in_nov-1)]
            last_laps_fam2_df=fam_df2[fam_df2[consts.LAP_COUNTER]>=frst_lap_fam2+laps_in_nov]
            if comp_to_last_only:
                last_laps_fam2_df=last_laps_fam2_df[(last_laps_fam2_df[consts.LAP_COUNTER]>=last_laps_fam2_df[consts.LAP_COUNTER].max()-laps_in_nov+1)]

            frst_laps_fam2_corrs=self.get_pairwise_correlation(params,frst_laps_fam2_df)
            last_laps_fam2_corrs=self.get_pairwise_correlation(params,last_laps_fam2_df)

            return frst_laps_fam2_corrs,last_laps_fam2_corrs
    
    

    





class DB_Analysis():
    def __init__(self):
        self.db_manager = DB_manager()

    def test_analysis_func(self, analysis_func, cage, mouse_name, seq, cell_num):
        exp = self.db_manager.get_single_experiment(cage, mouse_name, seq)
        cell = exp.cells[cell_num]
        result = analysis_func(cell)
        return result
    
    def test_analysis_func_on_exp(self, analysis_func, params, cage, mouse_name, seq):
        exp = self.db_manager.get_single_experiment(cage, mouse_name, seq)
        result = analysis_func(exp, params)
        return result

    def run_analysis_func(self, analysis_name, db_conditions, analysis_func, params, all=True, save=False):
        db = self.get_analysis_records(db_conditions, all)
        results = {}
        for row_num, record in db.iterrows():
            try:
                exp = self.db_manager.get_single_experiment(record[CAGE], record[MOUSE_NAME], record[SEQ])
            except:
                print(f'failed to load experiment for {record[CAGE]}_{record[MOUSE_NAME]}_{record[SEQ]}')
                continue
            for cell in exp.cells.values():
                # result = analysis_func(cell, params) # I deleted this for the mean time, im not sure where params should come into play. should every analysis function get params? 
                result = analysis_func(cell)
                results[cell] = (result)
        if save:
            self.save_results(results, analysis_name)
        return results
    
    def run_analysis_func_over_exps(self, analysis_name, db_conditions, analysis_func, params, all=True, save=False):
        db = self.get_analysis_records(db_conditions, all)
        results = {}
        for row_num, record in db.iterrows():
            try:
                exp = self.db_manager.get_single_experiment(record[CAGE], record[MOUSE_NAME], record[SEQ])
            except:
                print(f'failed to load experiment for {record[CAGE]}_{record[MOUSE_NAME]}_{record[SEQ]}')
                # break
            # result = analysis_func(exp, params) # I deleted this for the mean time, im not sure where params should come into play. should every analysis function get params? 
            result = analysis_func(exp)
            results[exp] = (result)
        if save:
            self.save_results(results, analysis_name)
        return results
    
        
    def get_analysis_records(self, db_conditions, all=True):
        db = self.db_manager.db
        if all:
            return db
        else:
            for param_name, param_values in db_conditions.items():
                if param_name[0]:  # if the first element is True, then we want to include the values
                    db = db.loc[db[param_name[1]].isin(param_values)] # param_name[1] is the name of the column
                else:
                    db = db.loc[~db[param_name[1]].isin(param_values)]
            return db

    def save_results(self, results, analysis_name):
        # svae results in a json format. the analysis name is the file name
        file_name = os.path.join(paths.ANALYSIS_RESULTS_DIR, analysis_name + '.json')
        with open(file_name, 'w') as fp:
            json.dump(results, fp, indent=4)

    
    def get_cell_names_by_cond(self, db_conditions):
        db = self.get_analysis_records(db_conditions, all)
    
        for param_name, param_values in db_conditions.items():
                if param_name[0]:  # if the first element is True, then we want to include the values
                    db = db.loc[db[param_name[1]].isin(param_values)] # param_name[1] is the name of the column
                else:
                    db = db.loc[~db[param_name[1]].isin(param_values)]
        cell_names=[]
        for row_num, record in db.iterrows():
            exp_name=record[CAGE]+'_'+record[MOUSE_NAME]+'_'+record[SEQ].split('_')[0]
            cells_as_string=record[GOOD_CELLS][1:-1]
            cell_nums=cells_as_string.split(',')
            for cell_num in cell_nums:
                cell_name=exp_name+'_'+cell_num.strip()
                cell_names.append(cell_name)
        return cell_names

            


    def does_test_result_exist_pkl(self, get_function):    #example usage: DB_Analysis().does_test_result_exist_pkl(Cell.get_spatial_info_Zscore)
        pkls_path = paths.ANALYSIS_RESULTS_DIR              # this function checks if a pkl file exists for a given get_function, under the assumption that the pkl file is named calc_function_name.pkl
        get_function_name = (str(get_function).split('.')[1]).split(' ')[0]
        calc_function_name = get_function_name.replace('get_', 'calc_')
        if os.path.exists(os.path.join(pkls_path, calc_function_name + '.pkl')):
            return True
        else:
            return False
        
    def load_data_from_pkl(self, get_function):
        pkls_path = paths.ANALYSIS_RESULTS_DIR
        get_function_name = (str(get_function).split('.')[1]).split(' ')[0]
        calc_function_name = get_function_name.replace('get_', 'calc_')
        with open(os.path.join(pkls_path, calc_function_name + '.pkl'), 'rb') as f:
            result_dict = pickle.load(f)
        return result_dict
    
    def get_data_from_pkl(self, get_function):
        if self.does_test_result_exist_pkl(get_function):
            return self.load_data_from_pkl(get_function)
        else:
            print('pkl file does not exist in folder')
            return None
    


    ################# results dictionary splitting functions #################

    def split_results_into_remapping_and_non_remapping(self, results_dict):
        # results is a dictionary of cell objects as keys and their results as values
        # return two dictionaries, one for remapping cells and one for non-remapping cells
        remapping_dict = {}
        non_remapping_dict = {}
        for cell, res in results_dict.items():
            if cell.metadata[REMAPPING]:
                remapping_dict[cell] = res
            else:
                non_remapping_dict[cell] = res
        return remapping_dict, non_remapping_dict
    
    def split_results_into_cell_types(self, results_dict):
        # results is a dictionary of cell objects as keys and their results as values
        # returns two dictionaries results by cell type
        pyr_dict = {}
        IN_dict = {}
        for cell, res in results_dict.items():
                if cell.metadata[CELL_TYPE] == 'Pyr':
                    pyr_dict[cell] = res
                else:
                    IN_dict[cell] = res
        return pyr_dict, IN_dict

    def split_results_within_a_cell(self,result_dict):
        # results is a dictionary of cell objects as keys and their results as values
        # return two lists, one for each index of the result
        res1 = []
        res2 = []
        for cell, res in result_dict.items():
                if res is not None: # if results are not nans, meaning the cell got a result 
                    res1.append(res[0])
                    res2.append(res[1])
        #turn nans into zeros
        return res1, res2 

    def split_results_into_cell_types_lists_no_remapping(self, results_dict):
        # results is a dictionary of cell objects as keys and their results as values
        # returns two lists results by cell type
        pyr_res = []
        IN_res = []
        for cell, res in results_dict.items():
            if ~cell.metadata[REMAPPING]:
                if cell.metadata[CELL_TYPE] == 'Pyr':
                    pyr_res.append(res[0])
                elif cell.metadata[CELL_TYPE] == 'IN':
                    IN_res.append(res[0])
                else:
                    print(f'cell {cell.name} has no cell type in form of Pyr or IN')
        return pyr_res, IN_res
    

    
    def split_results_into_cell_types_remapping(self, results_dict):
        # results is a dictionary of cell objects as keys and their results as values
        # returns four lists results by cell type and environment
        pyr_res_fam = []
        pyr_res_nov = []
        IN_res_fam = []
        IN_res_nov = []
        for cell, res in results_dict.items():
            if cell.metadata[REMAPPING]:
                if cell.metadata[CELL_TYPE] == 'Pyr':
                        pyr_res_fam.append(res[0])
                        pyr_res_nov.append(res[1])
                elif cell.metadata[CELL_TYPE] == 'IN':
                        IN_res_fam.append(res[0])
                        IN_res_nov.append(res[1])
                else:
                    print(f'cell {cell.name} has no cell type in form of Pyr or IN')
        #remove turn nan values into 0
        pyr_res_fam = np.nan_to_num(pyr_res_fam)
        pyr_res_nov = np.nan_to_num(pyr_res_nov)
        IN_res_fam = np.nan_to_num(IN_res_fam)
        IN_res_nov = np.nan_to_num(IN_res_nov)
        return pyr_res_fam, pyr_res_nov, IN_res_fam, IN_res_nov
    
    


    
    ############################# PLOTTING FUNCTIONS #########################################



    ################################## histogram plots ##################################

    def plot_result_histogram_by_cell_type(self, results, xlabel, bins=50,lims=None,alpha=0.8,pyr_color=PYR_COLOR,IN_color=IN_COLOR):
        # results is a dictionary of cell objects as keys and their results as values
        # analysis name is the name of the analysis, used for the title of the plot
        # lims is a list of two tuples, each tuple is the limits of the histogram for each cell type
        # plots the results in two different histograms, one for each cell type

        pyr_res, IN_res = self.split_results_into_cell_types_lists_no_remapping(results)

        plt.figure()
        if lims:
            plt.hist(pyr_res, label='Pyr', color=pyr_color, alpha=alpha, bins=bins, range=lims[0])
            plt.hist(IN_res, label='SST', color=IN_color, alpha=alpha, bins=bins, range=lims[1]);plt.ylabel('Count',fontsize=16);plt.legend(fontsize=14);plt.xlabel(xlabel,fontsize=16);plt.xticks(fontsize=16);plt.yticks(fontsize=16)
            plt.figure()
            plt.hist(pyr_res, label='Pyr', color=pyr_color, alpha=alpha, bins=bins, range=lims[0]);plt.ylabel('Count',fontsize=16);plt.legend(fontsize=14);plt.xlabel(xlabel,fontsize=16);plt.xticks(fontsize=16);plt.yticks(fontsize=16)
            plt.figure()
            plt.hist(IN_res, label='SST', color=IN_color, alpha=alpha, bins=bins, range=lims[1]);plt.ylabel('Count',fontsize=16);plt.legend(fontsize=14);plt.xlabel(xlabel,fontsize=16);plt.xticks(fontsize=16);plt.yticks(fontsize=16)
        else:
            plt.hist(pyr_res, label='Pyr', color=pyr_color, alpha=alpha, bins=bins)
            plt.hist(IN_res, label='SST', color=IN_color, alpha=alpha, bins=bins);plt.ylabel('Count',fontsize=16);plt.legend(fontsize=14);plt.xlabel(xlabel,fontsize=16);plt.xticks(fontsize=16);plt.yticks(fontsize=16)
            plt.figure()
            plt.hist(pyr_res, label='Pyr', color=pyr_color, alpha=alpha, bins=bins);plt.ylabel('Count',fontsize=16);plt.legend(fontsize=14);plt.xlabel(xlabel,fontsize=16);plt.xticks(fontsize=16);plt.yticks(fontsize=16)
            plt.figure()
            plt.hist(IN_res, label='SST', color=IN_color, alpha=alpha, bins=bins);plt.ylabel('Count',fontsize=16);plt.legend(fontsize=14);plt.xlabel(xlabel,fontsize=16);plt.xticks(fontsize=16);plt.yticks(fontsize=16)



    def plot_result_histogram_by_cell_type_remapping(self, results, analysis_name, bins=30,lims=None):
        pyr_res_fam, pyr_res_nov, IN_res_fam, IN_res_nov = self.split_results_into_cell_types_remapping(results)

        plt.figure()
        if lims:
            plt.hist(pyr_res_fam, label='Fam', color='r', alpha=0.4, bins=bins, range=lims[0])
            plt.hist(pyr_res_nov, label='Nov', color='m', alpha=0.4, bins=bins, range=lims[0])
        else:
            plt.hist(pyr_res_fam, label='Fam', color='r', alpha=0.4, bins=bins)
            plt.hist(pyr_res_nov, label='Nov', color='m', alpha=0.4, bins=bins)
        plt.title(analysis_name + ' ,Pyr')
        plt.legend()
        plt.show()
        
        plt.figure()
        if lims:
            plt.hist(IN_res_fam, label='Fam', color='b', alpha=0.4, bins=bins, range=lims[1])
            plt.hist(IN_res_nov, label='Nov', color='c', alpha=0.4, bins=bins, range=lims[1])
        else:
            plt.hist(IN_res_fam, label='Fam', color='b', alpha=0.4, bins=bins)
            plt.hist(IN_res_nov, label='Nov', color='c', alpha=0.4, bins=bins)
        plt.title(analysis_name + ' ,IN')
        plt.legend()
        plt.show()

        
    def plot_result_histogram(self, results, analysis_name, bins=50,limit=None):
        # plots all results in a histogram
        plt.figure()
        if limit:
            plt.hist(results.values(),color='k', bins=bins, range=limit)
        else:
            plt.hist(results.values(),color='k', bins=bins)
        plt.title(analysis_name)
        plt.show()

    ############################### end of histograms ########################################


    #################################### bar plots & violin plots #################################

    def plot_bar_plot_by_cell_type(self,results, analysis_name):
        # results is a dictionary of cell objects as keys and their results as values
        # plots the results in a bar plot, with error bars. pyrs in red, INs in blue
        pyr_res, IN_res = self.split_results_into_cell_types_lists(results)
        plt.figure()
        plt.bar(['Pyr', 'IN'], [np.mean(pyr_res), np.mean(IN_res)], yerr=[stats.sem(pyr_res), stats.sem(IN_res)], color=['r', 'b'], alpha=0.5)
        plt.title(analysis_name)
        plt.show()



    def plot_violin_plot_by_cell_type(self, result,ylabel):
        # results is a dictionary of cell objects as keys and their results as values
        # plots the results in a violin plot
        pyr_res, IN_res = self.split_results_into_cell_types_lists_no_remapping(result)
        pyr_res=np.array(pyr_res)
        pyr_res=pyr_res[pyr_res<500]
        sns.violinplot(data=[pyr_res, IN_res], palette=[PYR_COLOR, IN_COLOR],width=1), plt.xticks([0, 1], ['Pyr', 'SST'], fontsize=16), plt.ylabel(ylabel, fontsize=16)
        plt.yticks(fontsize=16)
        # preform wilcoxon rank sum test
        print(stats.ranksums(pyr_res, IN_res))
        # add significance stars
        x1, x2 = 0.1, 1
        y_max = max([max(pyr_res), max(IN_res)])
        y, h, col = y_max + y_max/50, y_max/50, 'k'
        plt.plot([x1, x1, x2-0.1, x2-0.1], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2-0.1)*.5, y+h, "***", ha='center', va='bottom', color=col, fontsize=20)

        plt.show()

    def plot_violin_plot_by_cell_type_remapping(self, result,ylabel):
        pyr_res_fam, pyr_res_nov, IN_res_fam, IN_res_nov=self.split_results_into_cell_types_remapping(result)
        #make two violin plots, one for pyrs and one for INs, comparing fam and nov
        plt.figure()
        plt.subplot(1,2,1)
        sns.violinplot(data=[pyr_res_fam, pyr_res_nov], palette=['r', 'm'],saturation=0.6), plt.xticks([0, 1], ['Fam', 'Nov'], fontsize=16), plt.ylabel(ylabel, fontsize=16)
        # preform wilcoxon rank sum test
        print('Pyr ',stats.ranksums(pyr_res_fam, pyr_res_nov))
        # preform paired t-test
        print('Pyr ', stats.ttest_rel(pyr_res_fam, pyr_res_nov))
        # add significance stars
        x1, x2 = 0.1, 0.9
        y_max = max([max(pyr_res_fam), max(pyr_res_nov)])
        y, h, col = y_max + y_max/50, y_max/50, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*.5, y+h, "n.s", ha='center', va='bottom', color=col, fontsize=20)
        plt.title('Pyr')
        plt.subplot(1,2,2)
        sns.violinplot(data=[IN_res_fam, IN_res_nov], palette=['b', 'c'],saturation=0.6), plt.xticks([0, 1], ['Fam', 'Nov'], fontsize=16)
        # preform wilcoxon rank sum test
        print('SST ',stats.ranksums(IN_res_fam, IN_res_nov))
        #perform paired t-test
        print('SST ',stats.ttest_rel(IN_res_fam, IN_res_nov))
        # add significance stars
        x1, x2 = 0.1, 0.9
        y_max = max([max(IN_res_fam), max(IN_res_nov)])
        y, h, col = y_max + y_max/50, y_max/50, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*.5, y+h, "***", ha='center', va='bottom', color=col, fontsize=20)
        plt.title('IN')
        plt.show()

    def plot_violin_plot_for_within_cell_comparison_of_two_vars(self,result,ylabel):
        #result is a dictionary of cell objects as keys and their results as values
        #plots the results in a violin plot
        first_var, second_var=self.split_results_within_a_cell(result)
        plt.figure()
        sns.violinplot(data=[first_var, second_var], palette=['b', 'c'],saturation=0.6), plt.xticks([0, 1], ['First Sixth', 'Last Sixth'], fontsize=16), plt.ylabel(ylabel, fontsize=16);
        # preform paired t-test
        print(stats.ttest_rel(first_var, second_var))
        # add significance stars
        x1, x2 = 0.1, 0.9
        y_max = max([max(first_var), max(second_var)])
        y, h, col = y_max + y_max/50, y_max/50, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*.5, y+h, "n.s", ha='center', va='bottom', color=col, fontsize=20)
        plt.show()


        ############################# end of bar plot & violin plot  #########################################


        ########################## per lap plots #########################################

    def plot_per_lap_results(self, result_dict, analysis_name, ylabel,laps):
        new_dict,new_dict_remap = self.leave_only_relevant_laps_and_arange_them(result_dict, 10)
        pyr_res, IN_res = self.split_results_into_cell_types(new_dict)
        pyr_res_remap, IN_res_remap = self.split_results_into_cell_types(new_dict_remap)
        pyr_stack=np.vstack(pyr_res.values())
        IN_stack=np.vstack(IN_res.values())
        pyr_stack_remap=np.vstack(pyr_res_remap.values())
        IN_stack_remap=np.vstack(IN_res_remap.values())    
        #plot the results
        plt.plot(np.mean(pyr_stack,0), 'r', label='Pyr');plt.legend();plt.title('No Remapping',fontsize=20);plt.ylabel('firing rate (Hz)',fontsize=16);plt.xlabel('Lap',fontsize=16)
        plt.ylim(0.8, 1.2)
        plt.show()
        plt.plot(np.mean(IN_stack,0), 'b', label='IN');plt.legend();plt.title('No Remapping',fontsize=20);plt.ylabel('firing rate (Hz)',fontsize=16);plt.xlabel('Lap',fontsize=16)
        plt.ylim(0.8, 1.2)
        plt.show()
        plt.plot(np.mean(pyr_stack_remap,0), 'r', label='Pyr');plt.legend();plt.title('With Remapping',fontsize=20);plt.ylabel('firing rate (Hz)',fontsize=16);plt.xlabel('Lap',fontsize=16)
        plt.ylim(0.8, 1.2)
        plt.show()
        plt.plot(np.mean(IN_stack_remap,0), 'b', label='IN');plt.legend();plt.title('With Remapping',fontsize=20);plt.ylabel('firing rate (Hz)',fontsize=16);plt.xlabel('Lap',fontsize=16)
        plt.ylim(0.8, 1.2)
        plt.show()

    def leave_only_relevant_laps_and_arange_them(self, results_dict, laps):
        new_dict = {}
        new_dict_remap = {}
        #laps is the number of laps to leave in the results
        # results is a dictionary of cell objects as keys and their results as values
        # returns a dictionary of cell objects as keys and their results as values with only the relevant laps
        for cell, res in results_dict.items():
            if len(res) == 1:
                if len(res[0]) >= laps*2:
                    middle_lap = int(len(res[0])/2)
                    fr_vector = res[0][middle_lap-laps:middle_lap+laps]
                    normalized_fr = fr_vector/np.mean(fr_vector)
                    new_dict[cell] = normalized_fr
            else:
                if (len(res[0]) >= laps) & (len(res[1]) >= laps):
                    fr_vector = np.concatenate((res[0][-laps:],res[1][0:laps]))
                    normalized_fr = fr_vector/np.mean(fr_vector)
                    new_dict_remap[cell] = normalized_fr
        return new_dict, new_dict_remap
    
        # result is a dictionary of cell objects as keys and their results as values
        # plots the results in a line plot
        


       # This class is for plotting the results of the analysis on many cells. It gets a dictionary of cell objects as keys and their results as values
# it should be able to plot the results in many different ways (bar plot, violin plot, histogram, etc.) it should also be able to split the results by different parameters (cell type, remap or no remap etc.)







class FigGenerator:
    def __init__(self, schema, db):
        self.schema = schema
        self.db = db
        self.CELL_NUMBER_PREFIX = "cell # "
        self.AX_TITLE_FRAMES = "Frames [#]"
        self.AX_TITLE_AU = "A.U"
        self.AX_TITLE_POSITION = 'Position'+'<br>'+'[cm]'
        self.AX_TITLE_FR = "Firing rate <br> [Hz]"
        self.CELL = "cell "
        self.SPIKES = " spikes"
        self.PYR = "Pyr"
        self.IN = "IN"
        self.CM_PER_BIN = " cm per bin"
        self.DELTA_F = "&#916;F/F"
        self.LOW_SEM = "low_sem"

    def get_longitundinal_metadata(self, experiment):
        long_metadata = self.db.loc[
            (self.db[CAGE] == experiment.metadata[CAGE]) & 
            (self.db[MOUSE_NAME] == experiment.metadata[MOUSE_NAME]) &
            (self.db[FOV] == experiment.metadata[FOV])
            ]
        long_metadata = long_metadata.sort_values(by=[EXPERIMENT_DATE], ascending=True)
        return long_metadata
            
    def get_longitudinal_data(self, experiment, cell_num):
        long_metadata = self.get_longitundinal_metadata(experiment)
        long_experiments = []
        for i, record in long_metadata.iterrows():
            if str(cell_num) in record[GOOD_CELLS]:
                db_record = {}
                for field, value in zip(self.schema, record):
                    if field == GOOD_CELLS:
                        value = str(value)
                    db_record[field] = value
                experiment = Experiment(self.schema, db_record)
                long_experiments.append(experiment)
        return long_experiments
    
    def get_generic_title(self, experiment):
        fig_title1 = experiment.metadata[CAGE] + " " + experiment.metadata[MOUSE_NAME] + " " + experiment.metadata[FOV]
        fig_title2 = "Pipeline seq: " + experiment.metadata[SEQ]
        return fig_title1 + "<br>" + fig_title2

    def get_axis_limit(self, trace, top_space, low_space):
        trace_ax_max = trace.max() + top_space * (trace.max() - trace.min())
        trace_ax_min = trace.min() - low_space * (trace.max() - trace.min())
        return trace_ax_min, trace_ax_max
    
    def get_range_edges_by_line_name(self, fig, line_name):
        """
        for generate button with two modes - relative scale
        and common scale to all subplots, dicts with edges 
        per each axis need to ne extracted
        """
        fig_min = float("inf")
        fig_max = float("-inf")
        relative_ranges = {}
        for line in fig.data:
            if line_name in line["name"]:
                ax_name = line["yaxis"]
                axis_key = 'yaxis.range'
                if len(ax_name) >= 1:
                    axis_key = 'yaxis{}.range'.format(ax_name[1:])
                ax_min = line["y"].min()
                ax_max = line["y"].max()
                relative_ranges[axis_key] = [ax_min, ax_max]
                if ax_max > fig_max:
                    fig_max = ax_max
                if ax_min < fig_min:
                    fig_min = ax_min
        common_range = {}
        for key in relative_ranges.keys():
            common_range[key] = [fig_min, fig_max]
        return relative_ranges, common_range

    def get_trace_color(self, experiment, fig_name=None):
        if experiment.metadata[CELL_TYPE]  == self.PYR:
            return 'red'
        if experiment.metadata[CELL_TYPE]  == self.IN:
            return 'darkblue'
        else:
            return 'red'
        
    def get_sem_color(self, experiment, fig_name=None):
        if experiment.metadata[CELL_TYPE]  == self.PYR:
            return 'indianred'
        if experiment.metadata[CELL_TYPE]  == self.IN:
            return 'lightblue'    
        else:
            return 'indianred'

    def get_spikes_color(self, experiment, fig_name=None):
        if experiment.metadata[CELL_TYPE]  == self.PYR:
            return 'blue'
        if experiment.metadata[CELL_TYPE]  == self.IN:
            return 'red'
        else:
            return 'blue'
        
    def get_subthreshold_color(self):
        trace_color = 'black'
        sem_color = 'gray'
        return trace_color, sem_color

    def get_bins_slider(self):
        return st.slider('Choose the number of bins:', 40, 200, 40, 20, )
    
    def get_cell_slider(self, experiment):
        if len(experiment.metadata[GOOD_CELLS]) == 1:
            return experiment.metadata[GOOD_CELLS][0]
        return st.select_slider('Choose cell number:',  options=experiment.metadata[GOOD_CELLS])

    def create_images_fig(self, fig_name, experiment):
        if fig_name is None:
            return self._create_image_fig_FOV(experiment)

    def _create_image_fig_FOV(self, experiment):
        volpy_data, slm_patterns, mean_image, _ = pipe_utils.get_pipline_results_data(
            experiment.metadata[CAGE], experiment.metadata[MOUSE_NAME], experiment.metadata[SEQ])
        images_fig = data_utils.display_FOV_images(mean_image, slm_patterns, volpy_data)
        return images_fig

    def create_fig(self, fig_name, experiment):
        if fig_name == CELLS_ACTIVITY:
            return self._create_fig_cells_activity(experiment)
        if fig_name == FR_AND_SUB:
            return self._create_fig_fr_and_sub(experiment)
        if fig_name == LAP_FR:
            return self._create_fig_lap_firing_rate(experiment)
        if fig_name == ACTIVITY_PER_LAP:
            return self._create_fig_activity_per_lap(experiment)
        if fig_name == LONGITUDIAL_ANALYSIS:
            return self._create_fig_longitudinal(experiment)
           
    def _create_fig_cells_activity(self, experiment):
        def _get_layout_cells_activity(experiment):
            rows_num = experiment.metadata[CELLS_NUM]
            row_heights = [25] * rows_num
            if experiment.behavior_flag:
                rows_num += 1
                row_heights += [5]
            fig = make_subplots \
                (
                rows = rows_num, cols = 1,
                vertical_spacing = 0.01, row_heights = row_heights, shared_xaxes=True, 
                row_titles = [self.CELL_NUMBER_PREFIX + str(i) for i in experiment.metadata[GOOD_CELLS]],
                y_title = self.DELTA_F , x_title = self.AX_TITLE_FRAMES
                )
            if experiment.behavior_flag:
                fig['layout']['yaxis'+str(experiment.metadata[CELLS_NUM]+1)]['title'] = self.AX_TITLE_POSITION
                fig['layout']['yaxis'+str(experiment.metadata[CELLS_NUM]+1)]['title'].update(font=dict(size=10))
                fig.update_yaxes(showticklabels=False, row=experiment.metadata[CELLS_NUM]+1, col=1)
            fig.update_layout(showlegend=True)
            fig.update_layout(title=self.get_generic_title(experiment), title_x=0.45, font=dict(size=18))
            return fig 
        
        def _add_data_cells_acivity(fig, experiment):
            df = experiment.preprocess_traces(experiment.raw_data)
            traces = experiment.get_traces(df)
            spikes_timming = experiment.get_spikes_timming(df)
            spike_heights = experiment.get_spikes_height(spikes_timming, traces)
            trace_color = self.get_trace_color(experiment)
            spike_color = self.get_spikes_color(experiment)
            for i, cell in enumerate(experiment.metadata[GOOD_CELLS]):
                trace = traces[i]
                spikes_times = spikes_timming[i]
                spikes_height = spike_heights[i]
                trace_ax_min, trace_ax_max = self.get_axis_limit(trace, 0.25, 0.25)
                fig.add_scatter(name=self.SPIKES, x=spikes_times, y=spikes_height, mode='markers', legendgroup =2, showlegend=i==0,
                                marker=dict(size=2.5, color=spike_color), visible='legendonly', row=i+1, col=1)
                fig.add_scatter(name=self.CELL +str(cell), x=np.arange(len(trace)), y=trace,
                                line=dict(color=trace_color, width=0.4), showlegend=False, row=i+1, col=1)
                fig.update_yaxes(row=i+1, col=1, range=[trace_ax_min, trace_ax_max])

            if experiment.behavior_flag:
                behave_col = experiment.get_position(df)
                lick_times = experiment.get_licks_timming(df)
                lick_value = 1.01 * behave_col[lick_times]
                fig.add_trace(go.Scatter(x=np.arange(len(behave_col)), y=behave_col, yaxis="y2",
                            showlegend=False, marker=dict(color='black')), row=experiment.metadata[CELLS_NUM]+1, col=1, )
                fig.add_scatter(name = " licks", x=lick_times, y=lick_value, mode='markers', marker=dict(
                    size=2.5, color=spike_color), visible='legendonly', row=experiment.metadata[CELLS_NUM]+1, col=1)
            return fig

        fig = _get_layout_cells_activity(experiment)
        fig = _add_data_cells_acivity(fig, experiment)
        return fig
    
    def _create_fig_fr_and_sub(self, experiment):
        bins_num = self.get_bins_slider()
        def _get_layout_fr_and_sub(experiment, bins_num):
            if experiment.metadata[REMAPPING]:
                cols = 2
                spacing_between_sub_plots = 0.1
                titels = ['<b>' + "FAMILIAR </b>" + "<br>Cell # " + str(experiment.metadata[GOOD_CELLS][0]),
                          '<b>' + "NOVEL </b>" + "<br>Cell # " + str(experiment.metadata[GOOD_CELLS][0])] + \
                        [item for sublist in [["Cell # " + str(i)]*2 for i in experiment.metadata[GOOD_CELLS][1:]] for item in sublist]
            else:
                cols = 1
                spacing_between_sub_plots = 0.05
                titels = [self.CELL_NUMBER_PREFIX + str(i) for i in experiment.metadata[GOOD_CELLS]]
            fig = make_subplots(
                rows = experiment.metadata[CELLS_NUM],
                cols = cols,
                shared_xaxes = True,
                vertical_spacing = spacing_between_sub_plots,
                specs=[[{"secondary_y": True}] * cols] * experiment.metadata[CELLS_NUM],
                subplot_titles=titels
                )

            cm_per_bin = experiment.calculate_bin_length(bins_num, experiment.data)
            fig_title = self.get_generic_title(experiment) 
            fig_title += '<br>Firing rate & subthreshold<br>' + str(cm_per_bin) + self.CM_PER_BIN
            fig.update_yaxes(title_text=self.AX_TITLE_FR, secondary_y=True, row=experiment.metadata[CELLS_NUM], col=cols,showgrid=False,zeroline=False)
            fig.update_yaxes(title_text=self.DELTA_F, secondary_y=False, row=experiment.metadata[CELLS_NUM], col=1,showgrid=False,zeroline=False)
            fig.update_xaxes(title_text=self.AX_TITLE_POSITION, row=experiment.metadata[CELLS_NUM],showgrid=False,zeroline=False)
            fig.update_xaxes(zeroline=False,showgrid=False)
            fig.update_yaxes(zeroline=False,showgrid=False)
            fig.update_layout(title=fig_title, title_x=0.45, title_y=0.975)
            fig.update_layout(paper_bgcolor='rgba(0,0,0,0)',  plot_bgcolor='rgba(0,0,0,0)')
            if experiment.metadata[CELLS_NUM] > 1:
                fig.update_layout(height=experiment.metadata[CELLS_NUM]*200, width=800)
            return fig

        def _add_data_fr_and_sub(fig, experiment, bins_num):
            rz_start_bin, rz_end_bin = experiment.get_reward_zone_bins(experiment.data, bins_num)
            if experiment.metadata[REMAPPING]:
                fam_df, nov_df = experiment.get_fam_and_novel_df()
                fam_df_fr = experiment.get_mean_firing_rate(bins_num, fam_df)
                nov_df_fr = experiment.get_mean_firing_rate(bins_num, nov_df)
                fam_df_sub = experiment.get_subthreshold_activity(bins_num, fam_df)
                nov_df_sub = experiment.get_subthreshold_activity(bins_num, nov_df)
                fam_df = fam_df_fr.merge(fam_df_sub, on=BIN, how='left')
                nov_df = nov_df_fr.merge(nov_df_sub, on=BIN, how='left')
                dfs = [fam_df, nov_df]
            else:
                df_fr = experiment.get_mean_firing_rate(bins_num)
                df_sub = experiment.get_subthreshold_activity(bins_num)
                df = df_fr.merge(df_sub, on=BIN, how='left')
                dfs = [df]
            trace_color = self.get_trace_color(experiment)
            sem_color = self.get_sem_color(experiment)
            subthresh_color, subthresh_sem_color = self.get_subthreshold_color()
            cm_per_bin = experiment.calculate_bin_length(bins_num, experiment.data)
            for col_num, df in enumerate(dfs, start=1):
                for row_num, cell_num in enumerate(experiment.metadata[GOOD_CELLS], start=1):
                    vis = True if row_num == 1 and col_num ==1 else False                    
                    fig.add_scatter(name='Firing rate', 
                                    x=df[BIN] * cm_per_bin, y=df[MEAN_FR_PREFIX + str(cell_num)], 
                                    row=row_num, col=col_num, 
                                    legendgroup=0, showlegend=vis, 
                                    line=dict(color=trace_color, width=4), secondary_y=True)
                    
                    fig.add_scatter(name=self.LOW_SEM,  
                                    x=df[BIN] * cm_per_bin,
                                    y=df[MEAN_FR_PREFIX + str(cell_num)]- df[SEM_FR_PREFIX + str(cell_num)], 
                                    mode='lines', line=dict(color=sem_color), secondary_y=True, 
                                    legendgroup=1, showlegend=False, fill=None, row=row_num, col=col_num, visible=True)
                    
                    fig.add_scatter(name='Firing rate SEM', 
                                    x=df[BIN] * cm_per_bin,
                                    y=df[MEAN_FR_PREFIX + str(cell_num)]+ df[SEM_FR_PREFIX + str(cell_num)], 
                                    mode='lines', line=dict(color=sem_color), secondary_y=True, legendgroup=1, 
                                    showlegend=vis, fill="tonexty", row=row_num, col=col_num, visible=True)
                    
                    fig.add_scatter(name="Subthreshold", 
                                    x=df[BIN] * cm_per_bin, y=smooth(df[SUB_ACTIVITY_PREFIX + str(cell_num)],3), 
                                    row=row_num, col=col_num, 
                                    line=dict(color=subthresh_color, width=3), 
                                    secondary_y=False, legendgroup=2, showlegend=vis, visible='legendonly')
                    
                    fig.add_scatter(name=self.LOW_SEM, 
                                      x=df[BIN] * cm_per_bin, 
                                      y=smooth(df[SUB_ACTIVITY_PREFIX + str(cell_num)] - df[SEM_SUB_ACTIVITY_PREFIX + str(cell_num)],3), 
                                      mode='lines', line=dict(color=subthresh_sem_color), secondary_y=False, 
                                      legendgroup=3, showlegend=False, fill=None, row=row_num, col=col_num, visible='legendonly')
                    
                    fig.add_scatter(name='Subthreshold SEM',
                                    x=df[BIN] * cm_per_bin,
                                    y=smooth(df[SUB_ACTIVITY_PREFIX + str(cell_num)] + df[SEM_SUB_ACTIVITY_PREFIX + str(cell_num)],3), 
                                    mode='lines', line=dict(color=subthresh_sem_color), secondary_y=False, 
                                    legendgroup=3, showlegend=vis, fill="tonexty", row=row_num, col=col_num, visible='legendonly')


                    fig.add_vrect(x0=rz_start_bin * cm_per_bin, x1=rz_end_bin * cm_per_bin, fillcolor="green", opacity=0.25,
                                line_width=0, row=row_num, col=col_num, annotation_font=dict(size=20, color="black"), 
                                annotation_text="<b>RZ</b>", annotation_position="top",)
            
            relative_ranges, common_range = self.get_range_edges_by_line_name(fig, 'Firing rate')
            scale_button = [dict(buttons=[
                    dict(label="relative_scale", method="relayout", args=[relative_ranges]),
                    dict(label="same_scale", method="relayout", args=[common_range])])]
            fig.update_layout(updatemenus=scale_button)

            return fig
        fig = _get_layout_fr_and_sub(experiment, bins_num)
        fig = _add_data_fr_and_sub(fig, experiment, bins_num)
        return fig 
    
    def _create_fig_lap_firing_rate(self, experiment):
        bins_num = self.get_bins_slider()
        def _get_layout_lap_firing_rate(experiment, bins_num):
            fig = make_subplots(
            rows=experiment.metadata[CELLS_NUM], cols=2,
            subplot_titles=tuple(
                [self.CELL_NUMBER_PREFIX + str(i) 
                for i in experiment.metadata[GOOD_CELLS]
                for _ in range(2)
                ]),  
            shared_xaxes=True, 
            vertical_spacing=0.1)
            fig.update_yaxes(title_text="Lap number", col=1)
            fig.update_xaxes(title_text="Position [cm]", row=experiment.metadata[CELLS_NUM])
            cm_per_bin = experiment.calculate_bin_length(bins_num, experiment.data)
            fig_title = self.get_generic_title(experiment) 
            fig_title += '<br>Firing rate per lap (' + str(cm_per_bin) + self.CM_PER_BIN + ')'
            fig.update_layout(title=fig_title, title_x=0.4)
            if experiment.metadata[CELLS_NUM] > 1:
                fig.update_layout(height=experiment.metadata[CELLS_NUM]*200, width=1000)
            else:
                fig.update_layout(height=400, width=1000)
            return fig
        def _add_data_lap_firing_rate(fig, experiment, bins_num):
            cm_per_bin = experiment.calculate_bin_length(bins_num, experiment.data)
            rz_start_bin, rz_end_bin = experiment.get_reward_zone_bins(experiment.data, bins_num)
            df = experiment.get_firing_rate_per_lap(bins_num)
            for i, cell_num in enumerate(experiment.metadata[GOOD_CELLS]):
                fig.add_trace(
                    go.Heatmap(
                    x=df[BIN] * cm_per_bin, 
                    y=df[consts.LAP_COUNTER],
                    z=df[FT_PREFIX + str(cell_num)], 
                    colorscale="amp", showscale=i == 0
                    ), 
                    row=i+1, col=2)
                
                fig.add_trace(
                    go.Heatmap(
                    x=df[BIN] * cm_per_bin, 
                    y=df[consts.LAP_COUNTER],
                    z=(df[FT_PREFIX + str(cell_num)] > 0).astype(int),
                    colorscale="amp", showscale=False
                    ),
                    row=i+1, col=1)

                fig.add_vrect(
                    x0=rz_start_bin * cm_per_bin,
                    x1=rz_end_bin * cm_per_bin, 
                    fillcolor="#c5d9ed", opacity=0.5,
                    line_width=0, row=i+1,  
                    annotation_font=dict(size=16, color="black"), 
                    annotation_text="<b>RZ</b>", 
                    annotation_position="top")

            return fig
        fig = _get_layout_lap_firing_rate(experiment, bins_num)
        fig = _add_data_lap_firing_rate(fig, experiment, bins_num)
        return fig 
    
    def _create_fig_activity_per_lap(self, experiment):
        cell_num = self.get_cell_slider(experiment)
        def _get_layout_activity_per_lap(experiment, cell_num):
            if experiment.metadata[REMAPPING]:
                cols_num = 2
                fam_df, nov_df = experiment.get_fam_and_novel_df()
                laps_num_fam = experiment.get_laps_number(fam_df)
                laps_num_nov = experiment.get_laps_number(nov_df)
                laps_num = max(laps_num_fam, laps_num_nov)
                width_scale = 2
                col_title = ["FAMILIAR", "NOVEL"]
            else:
                cols_num = 1
                laps_num = experiment.get_laps_number(experiment.data)
                width_scale = 1
                col_title = [""]
            fig = make_subplots(
                rows=laps_num, cols=cols_num, 
                shared_xaxes='all', shared_yaxes='all', 
                vertical_spacing=0.0, 
                row_width=[1/laps_num]*laps_num,
                row_titles=["lap #" + str(i+1) for i in range(laps_num)], 
                column_titles=col_title)
            if experiment.metadata[REMAPPING]:
                fig.update_xaxes(title_text=self.AX_TITLE_FRAMES, showticklabels=True, row=max(laps_num_fam, laps_num_nov))
            else:
                fig.update_xaxes(title_text=self.AX_TITLE_FRAMES, row=laps_num, col=1)
            fig.update_yaxes(showticklabels=False)
            fig.update_yaxes(showticklabels=True, row=laps_num//2, col=1)
            fig.update_yaxes(title_text=self.DELTA_F, row=laps_num//2, col=1)
            fig_title = self.get_generic_title(experiment) 
            fig_title += '<br>Activity per lap: Cell # ' + str(cell_num)
            fig.update_layout(title=fig_title, title_x=0.4)
            if experiment.metadata[CELLS_NUM] > 1:
                fig.update_layout(height=experiment.metadata[CELLS_NUM]*500, width=1000)
            else:
                fig.update_layout(height=1000, width=1000*width_scale)
            return fig
        
        def _add_data_activity_per_lap(fig, experiment, cell_num):
            if experiment.metadata[REMAPPING]:
                fam_df, nov_df = experiment.get_fam_and_novel_df()
                dfs = [fam_df, nov_df]
            else:
                dfs = [experiment.data]
            trace_col_prefix = data_utils.get_trace_col_prefix(dfs[0])
            trace_color = self.get_trace_color(experiment)
            for col_num, df in enumerate(dfs, start=1):
                for trace_num, (lap, lap_data) in enumerate(df.groupby([consts.LAP_COUNTER])):
                    trace = lap_data[trace_col_prefix + str(cell_num)]
                    lap_data = lap_data.reset_index()
                    lick_times = experiment.get_licks_timming(lap_data)
                    lick_value = trace.max() * np.ones(len(lick_times))
                    fig.add_scatter(name="env_"+str(col_num) + "_lap_" +str(lap), 
                                    x=np.arange(len(trace)), y=trace,
                                    line=dict(color=trace_color, width=0.4), 
                                    showlegend=False, 
                                    row=trace_num+1, col=col_num,
                                    xaxis='x'+str(col_num))

                    fig.add_scatter(name=" licks",
                                    x=lick_times, y=lick_value,
                                    mode='markers', 
                                    marker=dict(size=2.5, color='black'),
                                    visible='legendonly', 
                                    showlegend=(trace_num == 0) and (col_num == 1), 
                                    legendgroup=0, 
                                    row=trace_num+1, col=col_num, 
                                    xaxis='x'+str(col_num))

            return fig
        fig = _get_layout_activity_per_lap(experiment, cell_num)
        fig = _add_data_activity_per_lap(fig, experiment, cell_num)
        return fig 
    
    def _create_fig_longitudinal(self, experiment):       
        def _get_layout_longitudinal_fr_and_sub(experiment, cell_num):
            long_experiments = self.get_longitudinal_data(experiment, cell_num)
            sessions_num = len(long_experiments)

            remapping = False
            for exp in long_experiments:
                if exp.metadata[REMAPPING]:    
                    remapping = True

            if remapping:
                cols = 3
                column_widths = [0.2, 0.4, 0.4]
                column_titles=["Mean Image", "FAMILIAR", "NOVEL"]

            else:
                cols = 2
                column_widths=[0.3, 0.7]
                column_titles = ["Mean Image", "Firing Rate"]

            
            fig = make_subplots(
                rows=sessions_num, cols=cols, column_widths=column_widths, 
                shared_yaxes='columns', shared_xaxes='columns',
                specs=[[{"secondary_y": True}] * cols] * sessions_num,
                column_titles = column_titles,
                row_titles=["Week # " + str(i) for i in range(1, sessions_num+1)])

            fig.update_yaxes(title_text=self.AX_TITLE_FR, secondary_y=True, row=sessions_num, col=cols)
            fig.update_yaxes(title_text=self.DELTA_F, secondary_y=False, row=sessions_num, col=2)
            fig.update_xaxes(title_text=self.AX_TITLE_POSITION, row=sessions_num)

            fig.update_layout(coloraxis={'colorscale': 'gray'})
            fig.update_coloraxes(showscale=False)
            fig.update_yaxes(showticklabels=False,  col=1)
            fig.update_xaxes(showticklabels=False,  col=1)            
            fig_title = experiment.metadata[CAGE] + " " + \
            experiment.metadata[MOUSE_NAME] + " " + \
            experiment.metadata[FOV] + \
            " cell " + str(cell_num) + "<br>" + \
            str(sessions_num) + " Weeks of recordings"
            fig.update_layout(title=fig_title, title_x=0.5)
            for i in range(sessions_num):
                indent = 2
                if remapping:   
                    indent = 3
                fig.layout.annotations[i+indent].update(x=0,
                xshift=-65, textangle=0, font={'size': 12})
            # fig.layout.annotations[sessions_num+2].update(x=0.65)  # x title
            # fig.layout.annotations[sessions_num+3].update(x=0.38)  # y title
            # cancel share x axis on the first column
            fig.update_xaxes(matches=None, col=1)
            # cancel share x axis on the first column
            fig.update_yaxes(matches=None, col=1)
            fig.update_yaxes(autorange="reversed", col=1)
            return fig, long_experiments

        def _get_layout_longitudinal_activity(experiment, cell_num):
            long_experiments = self.get_longitudinal_data(experiment, cell_num)
            sessions_num = len(long_experiments)
            fig = make_subplots(
                rows=sessions_num, cols=2, column_widths=[0.3, 0.7], 
                shared_yaxes='columns', shared_xaxes='columns',
                y_title=self.DELTA_F, x_title="Time [ms]", 
                column_titles=["Mean Image", "Cell Activity"], 
                row_titles=["Week # " + str(i) for i in range(1, sessions_num+1)])
            fig.update_layout(coloraxis={'colorscale': 'gray'})
            fig.update_coloraxes(showscale=False)
            fig.update_yaxes(showticklabels=False,  col=1)
            fig.update_xaxes(showticklabels=False,  col=1)            
            fig_title = experiment.metadata[CAGE] + " " + \
            experiment.metadata[MOUSE_NAME] + " " + \
            experiment.metadata[FOV] + \
            " cell " + str(cell_num) + "<br>" + \
            str(sessions_num) + " Weeks of recordings"
            fig.update_layout(title=fig_title, title_x=0.5)
            for i in range(sessions_num):
                fig.layout.annotations[i+2].update(x=0,
                xshift=-65, textangle=0, font={'size': 12})
            fig.layout.annotations[sessions_num+2].update(x=0.65)  # x title
            fig.layout.annotations[sessions_num+3].update(x=0.38)  # y title
            # cancel share x axis on the first column
            fig.update_xaxes(matches=None, col=1)
            # cancel share x axis on the first column
            fig.update_yaxes(matches=None, col=1)
            fig.update_yaxes(autorange="reversed", col=1)
            return fig, long_experiments

        def _add_images_and_contours(fig, long_experiments, cell_num):
            for i, experiment in enumerate(long_experiments):
                
                _, slm_patterns, mean_image, _ = pipe_utils.get_pipline_results_data(
                    experiment.metadata[CAGE], experiment.metadata[MOUSE_NAME], experiment.metadata[SEQ])
                cell_contour = data_utils.get_cell_contour(slm_patterns, cell_num)

                fig_img = px.imshow(mean_image, origin='lower')
                fig.add_trace(fig_img.data[0], row=i+1, col=1)
                fig.add_trace(go.Scatter(
                    name="contour",
                    x=cell_contour[:, 0],
                    y=cell_contour[:, 1], 
                    line=dict(color="LightSkyBlue", width=2),
                    legendgroup=2, showlegend=i == 0,
                    ), row=i+1, col=1)
            return fig

        def _add_data_longitudinal_activity(fig, long_experiments, cell_num):            
            for i, experiment in enumerate(long_experiments):
                df = experiment.preprocess_traces(experiment.raw_data)
                traces = experiment.get_traces(df)
                spikes_timming = experiment.get_spikes_timming(df)
                spike_heights = experiment.get_spikes_height(spikes_timming, traces)
                trace_color = self.get_trace_color(experiment)
                spike_color = self.get_spikes_color(experiment)
                cell_idx = experiment.get_cell_idx(cell_num)
                trace_ax_min, trace_ax_max = self.get_axis_limit(traces[cell_idx], 0.25, 0.25)
                
                fig.add_scatter(name="Week # " +str(i+1), 
                                x=np.arange(len(traces[cell_idx])), 
                                y=traces[cell_idx],
                                line=dict(color=trace_color, width=0.4), 
                                showlegend=False, 
                                row=i+1, col=2)
                fig.add_scatter(name=self.SPIKES, 
                                x=spikes_timming[cell_idx], y=spike_heights[cell_idx], 
                                mode='markers', legendgroup =2, 
                                showlegend=i==0, visible='legendonly',
                                marker=dict(size=2.5, color=spike_color), 
                                row=i+1, col=2)
                fig.update_yaxes(row=i+1, col=2, range=[trace_ax_min, trace_ax_max])
            return fig
        
        def _add_data_longitudinal_fr_and_sub(fig, long_experiments, cell_num, bins_num):
            for row_num, experiment in enumerate(long_experiments, start=1):
                rz_start_bin, rz_end_bin = experiment.get_reward_zone_bins(experiment.data, bins_num)
                if experiment.metadata[REMAPPING]:
                    fam_df, nov_df = experiment.get_fam_and_novel_df()
                    fam_df_fr = experiment.get_mean_firing_rate(bins_num, fam_df)
                    nov_df_fr = experiment.get_mean_firing_rate(bins_num, nov_df)
                    fam_df_sub = experiment.get_subthreshold_activity(bins_num, fam_df)
                    nov_df_sub = experiment.get_subthreshold_activity(bins_num, nov_df)
                    fam_df = fam_df_fr.merge(fam_df_sub, on=BIN, how='left')
                    nov_df = nov_df_fr.merge(nov_df_sub, on=BIN, how='left')
                    dfs = [fam_df, nov_df]
                else:
                    df_fr = experiment.get_mean_firing_rate(bins_num)
                    df_sub = experiment.get_subthreshold_activity(bins_num)
                    df = df_fr.merge(df_sub, on=BIN, how='left')
                    dfs = [df]
                trace_color = self.get_trace_color(experiment)
                sem_color = self.get_sem_color(experiment)
                subthresh_color, subthresh_sem_color = self.get_subthreshold_color()
                cm_per_bin = experiment.calculate_bin_length(bins_num, experiment.data)
                for col_num, df in enumerate(dfs, start=2):
                    vis = True if row_num == 1 and col_num == 2 else False                    
                    fig.add_scatter(name='Firing rate', 
                                    x=df[BIN] * cm_per_bin, y=df[MEAN_FR_PREFIX + str(cell_num)], 
                                    row=row_num, col=col_num, 
                                    legendgroup=0, showlegend=vis, 
                                    line=dict(color=trace_color, width=4), secondary_y=True)
                    
                    fig.add_scatter(name=self.LOW_SEM,  
                                    x=df[BIN] * cm_per_bin,
                                    y=df[MEAN_FR_PREFIX + str(cell_num)]- df[SEM_FR_PREFIX + str(cell_num)], 
                                    mode='lines', line=dict(color=sem_color), secondary_y=True, 
                                    legendgroup=1, showlegend=False, fill=None, row=row_num, col=col_num, visible=True)
                    
                    fig.add_scatter(name='Firing rate SEM', 
                                    x=df[BIN] * cm_per_bin,
                                    y=df[MEAN_FR_PREFIX + str(cell_num)]+ df[SEM_FR_PREFIX + str(cell_num)], 
                                    mode='lines', line=dict(color=sem_color), secondary_y=True, legendgroup=1, 
                                    showlegend=vis, fill="tonexty", row=row_num, col=col_num, visible=True)
                    
                    fig.add_scatter(name="Subthreshold", 
                                    x=df[BIN] * cm_per_bin, y=df[SUB_ACTIVITY_PREFIX + str(cell_num)], 
                                    row=row_num, col=col_num, 
                                    line=dict(color=subthresh_color, width=3), 
                                    secondary_y=False, legendgroup=2, showlegend=vis, visible='legendonly')
                    
                    fig.add_scatter(name=self.LOW_SEM, 
                                        x=df[BIN] * cm_per_bin, 
                                        y=df[SUB_ACTIVITY_PREFIX + str(cell_num)] - df[SEM_SUB_ACTIVITY_PREFIX + str(cell_num)], 
                                        mode='lines', line=dict(color=subthresh_sem_color), secondary_y=False, 
                                        legendgroup=3, showlegend=False, fill=None, row=row_num, col=col_num, visible='legendonly')
                    
                    fig.add_scatter(name='Subthreshold SEM',
                                    x=df[BIN] * cm_per_bin,
                                    y=df[SUB_ACTIVITY_PREFIX + str(cell_num)] + df[SEM_SUB_ACTIVITY_PREFIX + str(cell_num)], 
                                    mode='lines', line=dict(color=subthresh_sem_color), secondary_y=False, 
                                    legendgroup=3, showlegend=vis, fill="tonexty", row=row_num, col=col_num, visible='legendonly')


                    fig.add_vrect(x0=rz_start_bin * cm_per_bin, x1=rz_end_bin * cm_per_bin, fillcolor="green", opacity=0.25,
                                line_width=0, row=row_num, col=col_num, annotation_font=dict(size=20, color="black"), 
                                annotation_text="<b>RZ</b>", annotation_position="top",)
            
                relative_ranges, common_range = self.get_range_edges_by_line_name(fig, 'Firing rate')
                scale_button = [dict(buttons=[
                        dict(label="relative_scale", method="relayout", args=[relative_ranges]),
                        dict(label="same_scale", method="relayout", args=[common_range])])]
                fig.update_layout(updatemenus=scale_button)
            return fig

        plot_type = st.radio("Choose a plot:", ('Activity', 'Firing rate'))
        cell_num = self.get_cell_slider(experiment)
        if plot_type == "Activity":
            fig, long_experiments = _get_layout_longitudinal_activity(experiment, cell_num)
            fig = _add_images_and_contours(fig, long_experiments, cell_num)
            fig = _add_data_longitudinal_activity(fig, long_experiments, cell_num)
        if plot_type == "Firing rate":
            bins_num = self.get_bins_slider()
            fig, long_experiments = _get_layout_longitudinal_fr_and_sub(experiment, cell_num)
            fig = _add_images_and_contours(fig, long_experiments, cell_num)
            fig = _add_data_longitudinal_fr_and_sub(fig, long_experiments, cell_num, bins_num)
        return fig 


    def _create_fig_template(self, experiment):
        bins_num = self.get_bins_slider()
        def _get_layout_template(experiment, bins_num):
            pass
            return fig
        def _add_data_template(fig, experiment, bins_num):
            pass
            return fig
        fig = _get_layout_template(experiment, bins_num)
        fig = _add_data_template(fig, experiment, bins_num)
        return fig

    

    

    
 
    

       


    
    
        


       

        


      
        

        


    



       


        

    


        

