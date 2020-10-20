# Ellen Hoffman fish data aggregation

## Pipeline 1: Best Recovery Drug
Goal: Multiple drugs were applied to mutant fish. Identify the one drug that will make the fish behave the most like the wild type.
Steps:
1. calculate zscore based on wild type, calculate_zscore.m/calculate_zscore_burst.m
2. average across all the fish for each drug and geno (HOM and WT), average_after_zscore.m
3. make clustergram, avgz_to_clustergram.m
4. PCA and euclidean distance.

Details:
The following two process lines work:

Individualdata ---> calculate_zscore ---> average_after_zscore ---> avgz_to_clustergram

splitmean ---> calculate_zscore_burst ---> average_after_zscore ---> avgz_to_clustergram

This one needs more work:

splitmean ---> calculate_zscore ---> average_after_zscore ---> avgz_to_clustergram

<table>
  <tr>
    <th>File name</th>
    <th>Description</th>
  </tr>
  <tr>
    <td>calculate_zscore.m</td>
    <td>Input type 1: Jeff's pipeline (usualy named 'individualData_xxxxxxxx.csv')<br>
      Output: calculate the z-score of experimental fish plate based on the mean and standard deviation of the wild type fish of the same day<br>
      Input type 2: The one off burst file generated by Jeff (named 'scn1lab_rw_split_means.csv')<br>
      Output: calculate the z-score of experimental fish plate based on the mean and standard deviation of the wild type fish of the same day, TO DO: removed burct and burur coloumns because there are too many infs</td>
  </tr>
  <tr>
    <td>calculate_zscore_burst.m</td>
    <td>Input: The one off burst file generated by Jeff (named 'scn1lab_rw_split_means.csv')<br>
      Output: calculate zscore of the burst variable (burct and burur), but use wild type from all days because each single given day the fish may not have any burst activity.</td>
  </tr>
  <tr>
    <td>average_after_zscore.m</td>
    <td>Input: calculated z score from the above 'calculate_zscore.m' and 'calculate_z_score_burst.m'<br>
      Trim the file so only HOM types are left<br>
      Average across fish, also aggregate across activities<br>
      Output1: mean_by_geno. The parameter/activities are averaged acrossed fish for each geno<br>
      If it is not the bust file, also export the aggregated zscore:  <br>
      Output2: _averaged. The parameter/activities are aggregated to rms and mean for each group, bout, activity, sleep, all</td>
  </tr>
  <tr>
    <td>avgz_to_clustergram.m</td>
    <td>Use the mean_by_geno to generate a clustergram.<br>
      Use customer color my_colormap. Replaced the underscores in the labels with space.</td>
  </tr>
  <tr>
    <td>run_PCA_euclidean.m/td>
    <td>Run PCA, make plots for publications, and calculate euclidean distances for all the drugs and dosages</td>
  </tr>
  <tr>
    <td>script_best_recovery_drug.m</td>
    <td>Combining the above steps to plots and calculate which drug is the best recovery drug.</td>
  </tr>
</table>

## Pipeline 2: Seizure analysis (pre-post) script_pre_post_analysis.m
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
    <td>plot the data by all the combinations of factors</td>
  </tr>
  <tr>
    <td>make_boxplots.m</td>
    <td>make boxplots for all the combinations of factors including geno, drug1, drug2, pre and post</td>
  </tr>
  <tr>
    <td>do_anova.m</td>
    <td>do 3-way anova on the difference score (geno, drug1, drug2), HOM and WT data (pre/post, drug1, drug2)</td>
  </tr>
  <tr>
    <td>script_pre_post_analysis.m</td>
    <td>a script that chains oher functions. First use  pre_post to create the output variable. Then use the output variable to make various plots and do anovas.</td>
  </tr>
</table>
