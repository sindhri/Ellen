# Ellen Hoffman fish data aggregation

The following two pipeline lines work:

Individualdata ---> calculate_zscore ---> average_after_zscore ---> avgz_to_clustergram

splitmean ---> calculate_zscore_burst ---> average_after_zscore ---> avgz_to_clustergram

This one needs more work:

splitmean ---> calculate_zscore ---> average_after_zscore ---> avgz_to_clustergram

## calculate_zscore.m
Input type 1: Jeff's pipeline (usualy named 'individualData_xxxxxxxx.csv')

Output: calculate the z-score of experimental fish plate based on the mean and standard deviation of the wild type fish of the same day

Input type 2: The one off burst file generated by Jeff (named 'scn1lab_rw_split_means.csv')

Output: calculate the z-score of experimental fish plate based on the mean and standard deviation of the wild type fish of the same day, TO DO: removed burct and burur coloumns because there are too many infs

## calculate_zscore_burst.m
Input: The one off burst file generated by Jeff (named 'scn1lab_rw_split_means.csv')

Output: calculate zscore of the burst variable (burct and burur), but use wild type from all days because each single given day the fish may not have any burst activity.

## average_after_zscore.m
Input: calculated z score from the above 'calculate_zscore.m' and 'calculate_z_score_burst.m'
Trim the file so only HOM types are left

Average across fish, also aggregate across activities

Output1: mean_by_geno. The parameter/activities are averaged acrossed fish for each geno
If it is not the bust file, also export the aggregated zscore:  
Output2: _averaged. The parameter/activities are aggregated to rms and mean for each group, bout, activity, sleep, all

## avgz_to_clustergram.m
Use the mean_by_geno to generate a clustergram.

Use customer color my_colormap. Replaced the underscores in the labels with space.

## Seizure analysis (pre-post) script_pre_post_analysis.m
* Summary: From Raw score for pre and post experiments to various plots, anovas, and an intermdeidate csv file
* Step 1: first make a folder that includes the following three files: pre, post, geno.
* Step 2: save the pre and post files into the excel format, with the xlsx extension. (Matlab is able to easily convert excel files into tables, but not from csv files)
* Step 3: run 'script_pre_post_analysis;' on the command line. You will select the three files in order: pre (the excel version), post (the excel version), geno. If it's an windows computer it would show the prompt for each type of file, but the prompts don't work in a Mac OS system. so just follow the order of pre-post-geno.
* Step 4: check the outputs which will be saved in the folder that you have created in the first place!

<table>
  <tr>
    <th>File name</th>
    <th>Description</th>
  </tr>
  <tr>
    <td>pre_post.m</td>
    <td>preprocess the pre and post data files and geno text file, generate an output that can feed into the following analyses</td>
  </tr>
  <tr>
    <td>plot_geno_by_time.m</td>
    <td>use the output from the pre_post and make geno_by_time plots, saved in a 'plots' folder</td>
  </tr>
  <tr>
    <td>plot_prepost.m</td>
    <td>plot the pre post data</td>
  </tr>
  <tr>
    <td>plot_bars.m</td>
    <td>plot the data by all the combinations of factors (to do: need to change to boxplots)</td>
  </tr>
  <tr>
    <td>do_anova.m</td>
    <td>do 3-way anova on the difference score: post-pre, (to do: split geno and do anova on both the pre and post data)</td>
  </tr>
  <tr>
    <td>script_pre_post_analysis.m</td>
    <td>a script that uses function pre_post and plot_geno_by_time, so that you can start from selecting the pre and post data file (in excel format), and geno text file, and the script would automatically generate geno_by_time plots.</td>
  </tr>
</table>
