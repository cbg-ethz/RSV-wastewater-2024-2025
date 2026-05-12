import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
from matplotlib.font_manager import FontProperties


# Assuming data is stored as a multi-line string or loaded from a CSV file
df_A = pd.read_csv('../../RSV/data_analysis/results/RSVA_2024_2025/relative_abundances/deconvolved.csv', sep='\t')
df_B = pd.read_csv('../../RSV/data_analysis/results/RSVB_2024_2025/relative_abundances/deconvolved.csv', sep='\t')
viral_loads_2024_2025 = pd.read_csv("/Users/arimaite/Documents/GitHub/combined_virus_sequencing/RSV/data_analysis/data/concentrations/rsv_ab_assay/viral_loads_imputated_rsvab.csv",
                                   usecols=['location', 'Date', 'RSV_A_imputed', 'RSV_B_imputed'])


def prep_df(df, location):
    df['location'] = df['location'].apply(lambda x: x.split(' ')[0])
    df['location'] = df['location'].replace({
        "Zürich (ZH)": "Zurich",
        "Genève (GE)": "Geneva"
    })
    df['location'] = df['location'].replace({
        "Zürich": "Zurich",
        "Genève": "Geneva"
    })

    df_filtered = df[df['location'] == location]

    df_filtered_new = pd.DataFrame()
     # make a new df with variants in different columns
    for variant in np.unique(df_filtered['variant']):
        new_col = df_filtered.loc[df_filtered['variant'] == variant, ['proportion', 'date']]
        df_filtered_new[variant] = new_col.reset_index()['proportion']

    time= pd.to_datetime(np.unique(df_filtered[['date']]))
    df_filtered_new.index = time
    # set a column with dates
    df_filtered_new['Date'] = df_filtered_new.index.values

    return df_filtered_new
def prep_viral_loads(viral_loads, location):
    # Extract viral loads for specific location
    data_rsv_viral_loads = viral_loads.loc[(viral_loads["location"] == location)].copy()
    data_rsv_viral_loads['Date'] = pd.to_datetime(data_rsv_viral_loads['Date'], format='%Y-%m-%d')
    # Sort viral loads by date
    data_rsv_viral_loads.sort_values(by='Date', inplace=True)

    return data_rsv_viral_loads


def plot_stacked_area_chart(df_A, df_B, location, filename, viral_loads):
    viral_loads_2024_2025_RSV_assay = viral_loads.loc[
        (viral_loads['RSV_A_imputed']).notna() & (viral_loads['RSV_B_imputed']).notna()].copy()


    viral_loads_2024_2025_RSV_assay['RSV_A_imputed'] = (
            viral_loads_2024_2025_RSV_assay['RSV_A_imputed']
    )
    viral_loads_2024_2025_RSV_assay['RSV_B_imputed'] = (
            viral_loads_2024_2025_RSV_assay['RSV_B_imputed']
    )

    df_filtered_new_RSVA = prep_df(df_A, location)
    df_filtered_new_RSVA.rename(columns = {
        'undetermined': 'undetermined_rsva'
    }, inplace=True)
    df_filtered_new_RSVB = prep_df(df_B, location)
    df_filtered_new_RSVB.rename(columns={
        'undetermined': 'undetermined_rsvb'
    }, inplace=True)
    #print(df_filtered_new_RSVA)
    #print(df_filtered_new_RSVB)
    rsva_variants = [
    "A.D.1",
    "A.D.1.4",
    "A.D.1.5",
    "A.D.3",
    "A.D.3.1",
    "A.D.5.2",
    "undetermined_rsva"]

    rsvb_variants = [
    "B.D.E.1",
    "B.D.E.1.1",
    "B.D.E.1.2",
    "B.D.E.1.8",
    "B.D.4.1.1",
    "B.D.E.7",
    "undetermined_rsvb"]

    df_merged = pd.merge(df_filtered_new_RSVA, df_filtered_new_RSVB, 'outer', on='Date')
    df_merged.sort_values(by='Date', inplace=True)
    df_merged['dropout_RSVA'] = np.nan
    df_merged['dropout_RSVB'] = np.nan

    df_merged['dropout_RSVA'] = np.where(df_merged[rsva_variants].isna().any(axis=1), 1, 0)
    df_merged['dropout_RSVB'] = np.where(df_merged[rsvb_variants].isna().any(axis=1), 1, 0)
    df_merged = df_merged.set_index('Date')
    #print(df_merged.shape)
    #print(df_merged.columns)
    data_rsv_viral_loads_RSV_prepared = prep_viral_loads(viral_loads_2024_2025_RSV_assay, location)
    #data_rsv_viral_loads_RSVB = prep_viral_loads(viral_loads_2024_2025_RSV_B, location)

    # extract only viral load values for the dates where deconvolution results exist
    #data_rsv_viral_loads_RSVA = data_rsv_viral_loads_RSVA.loc[data_rsv_viral_loads_RSVA['Date'].isin(df_merged.index)]
    #data_rsv_viral_loads_RSVB = data_rsv_viral_loads_RSVB.loc[data_rsv_viral_loads_RSVB['Date'].isin(df_merged.index)]

    # extract only deconvolution values for the dates where subtype-specific viral loads exist (i.e. subtyping was done)
    # TODO: interpolate concentration values where subtype-specific concentration was not measured -> DONE
    df_merged = df_merged.loc[df_merged.index.isin(data_rsv_viral_loads_RSV_prepared['Date'])]


# if the concentration proportion (RSV-AB assay) has been measured, then sample should be included even if there is no coverage!
    df_merged_old = df_merged

    df_merged = pd.merge(df_merged_old, data_rsv_viral_loads_RSV_prepared, 'outer', on='Date')



    df_merged = df_merged.iloc[:, :-4]

    df_merged["not_sequenced/dropout_RSVA"] = np.where(df_merged[rsva_variants].isna().any(axis=1), 1, 0)
    df_merged["not_sequenced/dropout_RSVB"] = np.where(df_merged[rsvb_variants].isna().any(axis=1), 1, 0)
    df_merged.sort_values(by='Date', inplace=True)
    df_merged.set_index('Date', inplace=True)

    df_merged.fillna(0, inplace=True)
    df_stacked_a = df_merged[rsva_variants+['not_sequenced/dropout_RSVA']].multiply(data_rsv_viral_loads_RSV_prepared.set_index('Date').reindex(df_merged.index)['RSV_A_imputed'], axis=0)


    df_stacked_b = df_merged[rsvb_variants+['not_sequenced/dropout_RSVB']].multiply(data_rsv_viral_loads_RSV_prepared.set_index('Date').reindex(df_merged.index)['RSV_B_imputed'], axis=0)
    df_stack_merged = pd.concat([df_stacked_a, df_stacked_b], axis=1)


    # Make the same indices in variants deconvolution df and viral loads df
    #data_rsv_viral_loads_RSVB.index = df_filtered_new_RSVB.index
    # Drop Date column
    #df_filtered_new_RSVB.drop(columns=['Date'], inplace=True)
    # Multiply proportions by total viral loads -> get abundance abundances for different variants (in the columns)
    #####df_filtered_new_stacked_area_RSVB = df_merged[rsvb_variants].mul(data_rsv_viral_loads_RSVB['total_RSVB_gc_ml'], axis=0)

    blues = [
        "navy",  # very dark blue
        "royalblue",  # strong medium blue
        "deepskyblue",  # bright, clear blue
        "steelblue",  # muted blue-grey
        "blueviolet",  # vivid purple-blue
        "mediumorchid",  # bright purple
        'black',
        'grey'
    ]

    greens = [
        "darkolivegreen",  # muted, darker green
        "seagreen",  # cooler green
        "limegreen",  # vivid bright green
        "mediumspringgreen",  # green with a mint tint
        "gold",  # strong yellow
        "khaki",
        "black",
        "grey"
    ]
    custom_colors = blues + greens

    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(10, 10))  # Adjust figure size for better readability
    #
    areas = df_stack_merged.T.values
    # #areas2 = df_filtered_new_stacked_area_RSVB.T.values
    # #x_vals = pd.to_datetime(df_filtered_new_stacked_area_RSVA.index) # df_filtered_new_stacked_area_RSVA
    #
    # # First group
    # #plot1 = ax.stackplot(x_vals, areas1, colors=[...], alpha=0.7)
    #
    # # Draw thick line between groups
    # #y_max = df_filtered_new_stacked_area_RSVA.sum(axis=1).values
    # #ax.plot(x_vals, y_max, color='black', linewidth=3)
    #
    # # Second group stacked on top
    # plot2 = ax.stackplot(df_merged.index, areas1,  alpha=0.7)
    #
    #
    # plt.show()
    #
    # Create the figure and axis
    # Plot the stackplot


    #print("labels:", df_stack_merged.columns)
    ax.stackplot(
        pd.to_datetime(df_stack_merged.index),  # x-axis values
        areas,                   # y-axis areas
        labels=(df_stack_merged.columns),  # Labels for each stack
        colors=custom_colors,    # Custom colors
        alpha=0.7                # Slight transparency for better visualization
    )
    y_max = df_stack_merged[rsva_variants+['not_sequenced/dropout_RSVA']].sum(axis=1).values
    ax.plot(pd.to_datetime(df_stack_merged.index), y_max, color='tomato', linewidth=2)
    ax.set_ylim(0, 3.1e07)

    time_monthly = pd.to_datetime(df_stack_merged.index).to_period('M').to_timestamp()

    # Drop duplicates to keep only one entry per month
    time_monthly = pd.to_datetime(np.unique(time_monthly))


    #ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=10, title="Variants")

    ax = plt.gca()
    ax.set_xticks(time_monthly)
    ax.set_xticklabels(
        (time_monthly).strftime("%b-%Y"),  # Ensure 'date' is a datetime type
        rotation=45,
        fontsize=25,
        ha='right')

    plt.xlabel('',fontsize=1)
    plt.ylabel('Load (gc/person/day)', fontsize=25)
    plt.title(f'{location}', fontsize=35)
    plt.legend(
        title='Variants',
        loc='upper center',  # center it horizontally
        bbox_to_anchor=(0.5, -0.4),  # position below the plot
        ncol=4,  # number of columns (adjust as needed)
        fontsize=12,
        title_fontsize=15,
        markerscale=0.8,
        handlelength=1.5,
        borderaxespad=0
    )
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # Force scientific notation
    ax.yaxis.get_offset_text().set_fontsize(14)  # increase size (default ~10)
    plt.xticks(rotation=45)
    ax.tick_params(axis='y', labelsize=25)

    plt.tight_layout()
    plt.savefig(filename)



plot_stacked_area_chart(df_A,
                        df_B,
                       location="Geneva",
                       filename="../../RSV/data_analysis/results/RSVB_2024_2025/relative_abundances/total_geneve_variants_over_time_stacked_area.pdf",
                       viral_loads=viral_loads_2024_2025)

plot_stacked_area_chart(df_A,
                        df_B,
                       location="Zurich",
                       filename="../../RSV/data_analysis/results/RSVB_2024_2025/relative_abundances/total_zurich_variants_over_time_stacked_area.pdf",
                       viral_loads=viral_loads_2024_2025)
plot_stacked_area_chart(df_A,
                        df_B,
                       location="Chur",
                       filename="../../RSV/data_analysis/results/RSVB_2024_2025/relative_abundances/total_rsvb_chur_variants_over_time_stacked_area.pdf",
                       viral_loads=viral_loads_2024_2025)
plot_stacked_area_chart(df_A,
                        df_B,
                       location="Basel",
                       filename="../../RSV/data_analysis/results/RSVB_2024_2025/relative_abundances/total_rsvb_basel_variants_over_time_stacked_area.pdf",
                       viral_loads=viral_loads_2024_2025)
plot_stacked_area_chart(df_A,
                        df_B,
                       location="Lugano",
                       filename="../../RSV/data_analysis/results/RSVB_2024_2025/relative_abundances/total_rsvb_lugano_variants_over_time_stacked_area.pdf",
                       viral_loads=viral_loads_2024_2025)
plot_stacked_area_chart(df_A,
                        df_B,
                       location="Laupen",
                       filename="../../RSV/data_analysis/results/RSVB_2024_2025/relative_abundances/total_rsvb_laupen_variants_over_time_stacked_area.pdf",
                       viral_loads=viral_loads_2024_2025)
