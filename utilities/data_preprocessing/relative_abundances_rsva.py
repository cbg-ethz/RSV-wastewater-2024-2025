import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
from matplotlib.font_manager import FontProperties


# Assuming data is stored as a multi-line string or loaded from a CSV file
df = pd.read_csv('../../RSV/data_analysis/results/RSVA_2024_2025/relative_abundances/deconvolved.csv', sep='\t')

# Function to create and save a plot for each location
def plot_by_location(df, location, filename):
    # Filter data by location

    df_filtered = df[df['location'] == location]

    df_filtered = df_filtered.sort_values(by=['date'])
    df_filtered['date'] = pd.to_datetime(df_filtered['date'])

    # Set the plot style
    sns.set_theme(style="whitegrid")
    #palette = sns.color_palette("tab20", n_colors=df_filtered['variant'].nunique())


    variant_colors = {
    "A.D.1": "navy",
    "A.D.1.4": "royalblue",
    "A.D.1.5": "deepskyblue",
    "A.D.3": "steelblue",
    "A.D.3.1": "blueviolet",
    "A.D.5.2": "mediumorchid",
    "undetermined": "black"}

    plt.figure(figsize=(20, 7))

    if location == "Zürich (ZH)":
        location = "Zurich (ZH)"
    elif location == "Genève (GE)":
        location = "Geneva (GE)"


    # Plot data for each variant
    for variant, group_data in df_filtered.groupby('variant'):

        plt.scatter(group_data['date'], group_data['proportion'], label=variant, color=variant_colors[variant])
        plt.plot(group_data['date'], group_data['proportion'], alpha=1, color=variant_colors[variant])
        plt.fill_between(group_data['date'], group_data['proportionLower'], group_data['proportionUpper'], alpha=0.1, color=variant_colors[variant])

        # # Add text for each point
        # for _, row in group_data.iterrows():
        #     value = row['proportion']
        #     lower = row['proportionLower']
        #     upper = row['proportionUpper']
        #     ci_text = f"{value:.2f}\n[{lower:.2f}, {upper:.2f}]"
        #     plt.text(
        #         row['date'],
        #         value + 0.02,  # small offset above the point
        #         ci_text,
        #         ha='center',
        #         va='bottom',
        #         fontsize=10,
        #         rotation=90  # optional, can set to 0 if you prefer horizontal
        #     )
        #



    # Customize the plot
    plt.xlabel('Date',fontsize=30)
    plt.ylabel('Proportion',fontsize=30)
    plt.title(f'{location}', fontsize=30)
    #plt.legend(title='Variants', bbox_to_anchor=(1.05, 1), loc='upper left',fontsize=20,title_fontsize=20)

    plt.gcf().subplots_adjust(bottom=0.2)  # Moves plot up to make space for x-axis label

    ax = plt.gca()  # Get the current Axes
    ax.tick_params(axis='y', labelsize=30)


    time= pd.to_datetime(np.unique(df_filtered[['date']]))
    plt.legend(title='Variants', bbox_to_anchor=(1.05, 1), loc='upper left',fontsize=30,title_fontsize=30, markerscale=2.5,
)
    # Take the first day of each month
    time_monthly = pd.to_datetime(time).to_period('M').to_timestamp()

    # Drop duplicates to keep only one entry per month
    time_monthly = pd.to_datetime(np.unique(time_monthly))
    #plt.xticks(rotation=45)

    ax = plt.gca()
    ax.set_xticks(time_monthly)
    ax.set_xticklabels(
        time_monthly.strftime("%b %Y"),  # you can format it like "2025.Apr"
        rotation=45,
        fontsize=30,
        ha='right')
    plt.tight_layout()


    # Save the plot
    plt.savefig(filename)


plot_by_location(df, "Genève (GE)", "../../RSV/data_analysis/results/RSVA_2024_2025/relative_abundances/rsva_geneve_variants_over_time.pdf")
plot_by_location(df, "Zürich (ZH)", "../../RSV/data_analysis/results/RSVA_2024_2025/relative_abundances/rsva_zurich_variants_over_time.pdf")
plot_by_location(df, "Laupen (BE)", "../../RSV/data_analysis/results/RSVA_2024_2025/relative_abundances/rsva_laupen_variants_over_time.pdf")
plot_by_location(df, "Chur (GR)", "../../RSV/data_analysis/results/RSVA_2024_2025/relative_abundances/rsva_chur_variants_over_time.pdf")
plot_by_location(df, "Lugano (TI)", "../../RSV/data_analysis/results/RSVA_2024_2025/relative_abundances/rsva_lugano_variants_over_time.pdf")
plot_by_location(df, "Basel (BS)", "../../RSV/data_analysis/results/RSVA_2024_2025/relative_abundances/rsva_basel_variants_over_time.pdf")

# def plot_stacked_area_chart(df, location, treatment_plant, filename, viral_loads):
#     df['location'] = df['location'].apply( lambda x: x.split(' ')[0])
#
#     df['location'] = df['location'].replace({
#         "Zürich (ZH)": "Zurich",
#         "Genève (GE)": "Geneva"
#     })
#     df['location'] = df['location'].replace({
#         "Zürich": "Zurich",
#         "Genève": "Geneva"
#     })
#     # Extract viral loads for specific location
#     data_rsv_viral_loads = viral_loads.loc[((viral_loads["Virus"] == "RSV-N")) & (viral_loads["location"] == treatment_plant)].copy()
#     data_rsv_viral_loads['Date'] = pd.to_datetime(data_rsv_viral_loads['Date'], format='%Y-%m-%d')
#     # Sort viral loads by date
#     data_rsv_viral_loads.sort_values(by='Date', inplace=True)
#     # extract relative abundances from specific location
#     df_filtered = df[df['location'] == location]
#
#
#     df_filtered_new = pd.DataFrame()
#      # make a new df with variants in different columns
#     for variant in np.unique(df_filtered['variant']):
#         new_col = df_filtered.loc[df_filtered['variant'] == variant, ['proportion','date']]
#         df_filtered_new[variant] = new_col.reset_index()['proportion']
#
#     time= pd.to_datetime(np.unique(df_filtered[['date']]))
#     df_filtered_new.index = time
#     # set a column with dates
#     df_filtered_new['Date'] = df_filtered_new.index.values
#
#     # extract only viral load values for the dates where deconvolution results exist
#     data_rsv_viral_loads = data_rsv_viral_loads.loc[data_rsv_viral_loads['Date'].isin(df_filtered_new.index)]
#     # extract only deconvolution values for the dates where viral loads exist
#     df_filtered_new = df_filtered_new.loc[df_filtered_new['Date'].isin(data_rsv_viral_loads['Date'])]
#
#     # Make the same indices in variants deconvolution df and viral loads df
#     data_rsv_viral_loads.index = df_filtered_new.index
#     # Drop Date column
#     df_filtered_new.drop(columns=['Date'], inplace=True)
#     # Multiply proportions by total viral loads -> get abundance abundances for different variants (in the columns)
#     df_filtered_new_stacked_area = df_filtered_new.mul(data_rsv_viral_loads['total_RSVA_gc_ml'], axis=0)
#     print(df_filtered_new_stacked_area.head())
# # Define custom colors for the plot
#     custom_colors = [
#         "navy",  # very dark blue
#         "royalblue",  # strong medium blue
#         "deepskyblue",  # bright, clear blue
#         "steelblue",  # muted blue-grey
#         "blueviolet",  # vivid purple-blue
#         "mediumorchid",  # bright purple
#         "black"
#
#     ]
#
#
#     sns.set_theme(style="whitegrid")
#
#     areas = df_filtered_new_stacked_area.T.values
#
#     print("areas shape:", areas.shape)
#
#     # Create the figure and axis
#     fig, ax = plt.subplots(figsize=(20, 10))  # Adjust figure size for better readability
#     # Plot the stackplot
#
#     ax.stackplot(
#         pd.to_datetime(df_filtered_new_stacked_area.index),  # x-axis values
#         areas,                   # y-axis areas
#         labels=(df_filtered_new_stacked_area.columns),  # Labels for each stack
#         colors=custom_colors,    # Custom colors
#         alpha=0.7                # Slight transparency for better visualization
#     )
#     ax.set_ylim(0, 2*10**7)
#
#     time_monthly = pd.to_datetime(time).to_period('M').to_timestamp()
#     print(time_monthly)
#
#     # Drop duplicates to keep only one entry per month
#     time_monthly = pd.to_datetime(np.unique(time_monthly))
#
#
#     ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=30, title="Variants")
#
#     ax = plt.gca()
#     ax.set_xticks(time_monthly)
#     ax.set_xticklabels(
#         (time_monthly).strftime("%Y.%b.%d"),  # Ensure 'date' is a datetime type
#         rotation=45,
#         fontsize=30,
#         ha='right')
#
#     plt.xlabel('Date',fontsize=20)
#     plt.ylabel('Flow-normalized viral load \n (gc/person/day)',fontsize=20)
#     plt.title(f'2024-2025 {location}', fontsize=25)
#     plt.legend(title='Variants', bbox_to_anchor=(1.05, 1), loc='upper left',fontsize=20,title_fontsize=20)
#     ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
#     ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # Force scientific notation
#     plt.xticks(rotation=45)
#     ax.tick_params(axis='y', labelsize=30)
#
#     plt.tight_layout()
#     plt.savefig(filename)
# #
# # plot_stacked_area_chart(df=df,
# #                        location="Geneva",
# #                        treatment_plant ="Geneva",
# #                        filename="../../RSV/data_analysis/results/RSVA_2024_2025/relative_abundances/rsva_geneve_variants_over_time_stacked_area.pdf",
# #                        viral_loads=viral_loads_2024_2025)
# #
# # plot_stacked_area_chart(df=df,
# #                        location="Zurich",
# #                        treatment_plant ="Zurich",
# #                        filename="../../RSV/data_analysis/results/RSVA_2024_2025/relative_abundances/rsva_zurich_variants_over_time_stacked_area.pdf",
# #                        viral_loads=viral_loads_2024_2025)
# # plot_stacked_area_chart(df=df,
# #                        location="Chur",
# #                        treatment_plant ="Chur",
# #                        filename="../../RSV/data_analysis/results/RSVA_2024_2025/relative_abundances/rsva_chur_variants_over_time_stacked_area.pdf",
# #                        viral_loads=viral_loads_2024_2025)
# # plot_stacked_area_chart(df=df,
# #                        location="Basel",
# #                        treatment_plant ="Basel",
# #                        filename="../../RSV/data_analysis/results/RSVA_2024_2025/relative_abundances/rsva_basel_variants_over_time_stacked_area.pdf",
# #                        viral_loads=viral_loads_2024_2025)
# # plot_stacked_area_chart(df=df,
# #                        location="Lugano",
# #                        treatment_plant ="Lugano",
# #                        filename="../../RSV/data_analysis/results/RSVA_2024_2025/relative_abundances/rsva_lugano_variants_over_time_stacked_area.pdf",
# #                        viral_loads=viral_loads_2024_2025)
# # plot_stacked_area_chart(df=df,
# #                        location="Laupen",
# #                        treatment_plant ="Laupen",
# #                        filename="../../RSV/data_analysis/results/RSVA_2024_2025/relative_abundances/rsva_laupen_variants_over_time_stacked_area.pdf",
# #                        viral_loads=viral_loads_2024_2025)