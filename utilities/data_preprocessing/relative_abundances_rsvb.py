import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.ticker as ticker

# Assuming data is stored as a multi-line string or loaded from a CSV file
df = pd.read_csv('../../RSV/data_analysis/results/RSVB_2024_2025/relative_abundances/deconvolved.csv', sep='\t')

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
    "B.D.E.1": "darkolivegreen",
    "B.D.E.1.1": "seagreen",
    "B.D.E.1.2": "limegreen",
    "B.D.E.1.8": "mediumspringgreen",
    "B.D.4.1.1": "gold",
    "B.D.E.7": "khaki",
    "undetermined": "black"}

    if location == "Zürich (ZH)":
        location = "Zurich (ZH)"
    elif location == "Genève (GE)":
        location = "Geneva (GE)"


    #variant_colors = dict(zip(df_filtered['variant'].unique(), palette))
    plt.figure(figsize=(20, 7))

    for variant, group_data in df_filtered.groupby('variant'):

        plt.scatter(group_data['date'], group_data['proportion'], label=variant, color=variant_colors[variant])
        plt.plot(group_data['date'], group_data['proportion'], alpha=1, color=variant_colors[variant])
        plt.fill_between(group_data['date'], group_data['proportionLower'], group_data['proportionUpper'], alpha=0.1, color=variant_colors[variant])

        # Add text for each point
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


    # Customize the plot
    plt.xlabel('Date',fontsize=30)
    plt.ylabel('Proportion',fontsize=30)
    plt.title(f'{location}', fontsize=30)

    plt.gcf().subplots_adjust(bottom=0.2)  # Moves plot up to make space for x-axis label

    ax = plt.gca()  # Get the current Axes
    ax.tick_params(axis='y', labelsize=30)



    time= pd.to_datetime(np.unique(df_filtered[['date']]))
    plt.legend(title='Variants', bbox_to_anchor=(1.05, 1), loc='upper left',fontsize=30,title_fontsize=30, markerscale=2.5)
    # Take the first day of each month
    time_monthly = pd.to_datetime(time).to_period('M').to_timestamp()

    # Drop duplicates to keep only one entry per month
    time_monthly = pd.to_datetime(np.unique(time_monthly))
    #plt.xticks(rotation=45)

    ax = plt.gca()
    ax.set_xticks(time_monthly)
    ax.set_xticklabels(
        time_monthly.strftime("%b %Y"),
        rotation=45,
        fontsize=30,
        ha='right')
    plt.tight_layout()

    # plt.show()
    plt.savefig(filename)


plot_by_location(df, "Genève (GE)", "../../RSV/data_analysis/results/RSVB_2024_2025/relative_abundances/rsvb_geneve_variants_over_time.pdf")
plot_by_location(df, "Zürich (ZH)", "../../RSV/data_analysis/results/RSVB_2024_2025/relative_abundances/rsvb_zurich_variants_over_time.pdf")
plot_by_location(df, "Laupen (BE)", "../../RSV/data_analysis/results/RSVB_2024_2025/relative_abundances/rsvb_laupen_variants_over_time.pdf")
plot_by_location(df, "Chur (GR)", "../../RSV/data_analysis/results/RSVB_2024_2025/relative_abundances/rsvb_chur_variants_over_time.pdf")
plot_by_location(df, "Lugano (TI)", "../../RSV/data_analysis/results/RSVB_2024_2025/relative_abundances/rsvb_lugano_variants_over_time.pdf")
plot_by_location(df, "Basel (BS)", "../../RSV/data_analysis/results/RSVB_2024_2025/relative_abundances/rsvb_basel_variants_over_time.pdf")

