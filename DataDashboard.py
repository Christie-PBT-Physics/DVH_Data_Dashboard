import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
import os
import re
import collections


@st.cache_data
def read_csv(location):
    arr = []
    for file in os.listdir(location):
        match = re.search(r'(\d*)_All_DVH_Data.csv',file)
        if (match):
            PID = match[1]
            full_path = os.path.join(directory, file)
            try:
                with open(full_path, 'r') as file:
                    organ_at_risk = ''
                    max_dose = 0
                    mean_dose = 0
                    meadian_dose = 0
                    min_dose = 99999
                    organ_volume = 0
                    found_oar = False
                    found_metrics = False
                    for line in file:
                        if not found_oar:
                            if len(line.split(r',')) == 1:
                                organ_at_risk = line.strip('\n')
                                found_oar = True
                        else:
                            if not found_metrics:
                                if line == 'Max,Mean,Median,Min,Volume\n':
                                    found_metrics = True
                            else:
                                data_row = line.split(',')
                                data_row = [float(value.strip(' cGy')) for value in data_row]
                                max_dose, mean_dose, meadian_dose, min_dose, organ_volume = data_row
                                found_oar = False
                                found_metrics = False
                                arr.append([PID, organ_at_risk, max_dose, mean_dose, meadian_dose, min_dose, organ_volume])
            except:
                print(f'Could not find {path}')
                    
    df = pd.DataFrame(arr, columns=('Patient', 'OAR', 'Max', 'Mean' ,'Median' , 'Min', 'Volume'))
    number = len(df['Patient'].unique())
    return df, number


def plot_single_violin(all, patient, metric, oar):
    plot = plt.figure(figsize=(10,1), dpi=300)
    sns.violinplot(data=all[all['OAR']==oar],
                    y='OAR',
                    x=metric,
                    cut=0,
                    bw_adjust=.5,
                    linewidth=0,
                    inner="stick",
                    legend=False,
                    inner_kws=dict(linewidth=.8, color=".8"))
    sns.stripplot(data=patient[patient['OAR']==oar], y='OAR', x=metric, color='red')
    plt.yticks([])
    plt.ylabel('')
    plt.xlabel('')
    plt.xticks(fontsize=8)
    plt.box(False)
    plt.tight_layout()
    plt.grid(True, linestyle=':', linewidth=.5)
    return plot


def plot_single_violin_plotly(all_data, patient, metric, oar):
    specific_oar_df = all_data[all_data['OAR']==oar]
    violin = go.Violin(x=specific_oar_df[metric],
                        points=False,
                        # jitter=1,
                        # pointpos=0,
                        marker=dict(color='#4575b4'),
                        spanmode='hard',
                        bandwidth=1,
                        name='violin',
                        hoverinfo='none',
                        zorder=2)
    scatter = go.Scatter(x=patient[patient['OAR']==oar][metric],
                            y=['violin'],
                            marker=dict(color='#a50026'),
                            zorder=3)
    points_df = swarm_df(all_data[all_data['OAR']==oar][metric])
    points = go.Scatter(x=points_df['x'],
                            y=points_df['y'],
                            marker=dict(color='#abd9e9', symbol='circle-open'),
                            mode='markers',
                            yaxis='y2',
                            zorder=1)

    fig = go.Figure(data=[scatter,points,violin])
    
    fig.update_layout(yaxis=dict(visible=False),
                        yaxis2=dict(visible=False, overlaying='y', side='right'),
                        xaxis=dict(showgrid=True),
                        margin=dict(l=20, r=20, t=20, b=20),
                        height=150,
                        showlegend=False)
    return fig


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# https://plotly.com/python/line-and-scatter/ (Swarm (or Beeswarm) Plots)
def negative_1_if_count_is_odd(count):
    # if this is an odd numbered entry in its bin, make its y coordinate negative
    # the y coordinate of the first entry is 0, so entries 3, 5, and 7 get
    # negative y coordinates
    if count % 2 == 1:
        return -1
    else:
        return 1

def swarm(
    X_series,
    fig_title,
    point_size=2,
    fig_width=800,
    gap_multiplier=1.2,
    bin_fraction=0.95,  # slightly undersizes the bins to avoid collisions
):
    # sorting will align columns in attractive c-shaped arcs rather than having
    # columns that vary unpredictably in the x-dimension.
    # We also exploit the fact that sorting means we see bins sequentially when
    # we add collision prevention offsets.
    X_series = X_series.copy().sort_values()

    # we need to reason in terms of the marker size that is measured in px
    # so we need to think about each x-coordinate as being a fraction of the way from the
    # minimum X value to the maximum X value
    min_x = min(X_series)
    max_x = max(X_series)

    list_of_rows = []
    # we will count the number of points in each "bin" / vertical strip of the graph
    # to be able to assign a y-coordinate that avoids overlapping
    bin_counter = collections.Counter()

    for x_val in X_series:
        # assign this x_value to bin number
        # each bin is a vertical strip slightly narrower than one marker
        bin = (((fig_width*bin_fraction*(x_val-min_x))/(max_x-min_x)) // point_size)

        # update the count of dots in that strip
        bin_counter.update([bin])

        # remember the "y-slot" which tells us the number of points in this bin and is sufficient to compute the y coordinate unless there's a collision with the point to its left
        list_of_rows.append(
            {"x": x_val, "y_slot": bin_counter[bin], "bin": bin})

    # iterate through the points and "offset" any that are colliding with a
    # point to their left apply the offsets to all subsequent points in the same bin.
    # this arranges points in an attractive swarm c-curve where the points
    # toward the edges are (weakly) further right.
    bin = 0
    offset = 0
    for row in list_of_rows:
        if bin != row["bin"]:
            # we have moved to a new bin, so we need to reset the offset
            bin = row["bin"]
            offset = 0
        # see if we need to "look left" to avoid a possible collision
        for other_row in list_of_rows:
            if (other_row["bin"] == bin-1):
                # "bubble" the entry up until we find a slot that avoids a collision
                while ((other_row["y_slot"] == row["y_slot"]+offset) and (((fig_width*(row["x"]-other_row["x"]))/(max_x-min_x) // point_size) < 1)):
                    offset += 1
                    # update the bin count so we know whether the number of
                    # *used* slots is even or odd
                    bin_counter.update([bin])

        row["y_slot"] += offset
        # The collision free y coordinate gives the items in a vertical bin
        # y-coordinates to evenly spread their locations above and below the
        # y-axis (we'll make a correction below to deal with even numbers of
        # entries).  For now, we'll assign 0, 1, -1, 2, -2, 3, -3 ... and so on.
        # We scale this by the point_size*gap_multiplier to get a y coordinate
        # in px.
        row["y"] = (row["y_slot"]//2) * \
            negative_1_if_count_is_odd(row["y_slot"])*point_size*gap_multiplier

    # if the number of points is even, move y-coordinates down to put an equal
    # number of entries above and below the axis
    for row in list_of_rows:
        if bin_counter[row["bin"]] % 2 == 0:
            row["y"] -= point_size*gap_multiplier/2

    df = pd.DataFrame(list_of_rows)
    # One way to make this code more flexible to e.g. handle multiple categories
    # would be to return a list of "swarmified" y coordinates here and then plot
    # outside the function.
    # That generalization would let you "swarmify" y coordinates for each
    # category and add category specific offsets to put the each category in its
    # own row

    fig = px.scatter(
        df,
        x="x",
        y="y",
        title=fig_title,
    )
    # we want to suppress the y coordinate in the hover value because the
    # y-coordinate is irrelevant/misleading
    fig.update_traces(
        marker_size=point_size,
        # suppress the y coordinate because the y-coordinate is irrelevant
        hovertemplate="<b>value</b>: %{x}",
    )
    # we have to set the width and height because we aim to avoid icon collisions
    # and we specify the icon size in the same units as the width and height
    fig.update_layout(width=fig_width, height=(
        point_size*max(bin_counter.values())+200))
    fig.update_yaxes(
        showticklabels=False,  # Turn off y-axis labels
        ticks='',               # Remove the ticks
        title=""
    )
    return fig

def swarm_df(
    X_series,
    point_size=0.5,
    fig_width=800,
    gap_multiplier=1.1,
    bin_fraction=0.95,  # slightly undersizes the bins to avoid collisions
):
    # sorting will align columns in attractive c-shaped arcs rather than having
    # columns that vary unpredictably in the x-dimension.
    # We also exploit the fact that sorting means we see bins sequentially when
    # we add collision prevention offsets.
    X_series = X_series.copy().sort_values()

    # we need to reason in terms of the marker size that is measured in px
    # so we need to think about each x-coordinate as being a fraction of the way from the
    # minimum X value to the maximum X value
    min_x = min(X_series)
    max_x = max(X_series)

    list_of_rows = []
    # we will count the number of points in each "bin" / vertical strip of the graph
    # to be able to assign a y-coordinate that avoids overlapping
    bin_counter = collections.Counter()

    for x_val in X_series:
        # assign this x_value to bin number
        # each bin is a vertical strip slightly narrower than one marker
        bin = (((fig_width*bin_fraction*(x_val-min_x))/(max_x-min_x)) // point_size)

        # update the count of dots in that strip
        bin_counter.update([bin])

        # remember the "y-slot" which tells us the number of points in this bin and is sufficient to compute the y coordinate unless there's a collision with the point to its left
        list_of_rows.append(
            {"x": x_val, "y_slot": bin_counter[bin], "bin": bin})

    # iterate through the points and "offset" any that are colliding with a
    # point to their left apply the offsets to all subsequent points in the same bin.
    # this arranges points in an attractive swarm c-curve where the points
    # toward the edges are (weakly) further right.
    bin = 0
    offset = 0
    for row in list_of_rows:
        if bin != row["bin"]:
            # we have moved to a new bin, so we need to reset the offset
            bin = row["bin"]
            offset = 0
        # see if we need to "look left" to avoid a possible collision
        for other_row in list_of_rows:
            if (other_row["bin"] == bin-1):
                # "bubble" the entry up until we find a slot that avoids a collision
                while ((other_row["y_slot"] == row["y_slot"]+offset) and (((fig_width*(row["x"]-other_row["x"]))/(max_x-min_x) // point_size) < 1)):
                    offset += 1
                    # update the bin count so we know whether the number of
                    # *used* slots is even or odd
                    bin_counter.update([bin])

        row["y_slot"] += offset
        # The collision free y coordinate gives the items in a vertical bin
        # y-coordinates to evenly spread their locations above and below the
        # y-axis (we'll make a correction below to deal with even numbers of
        # entries).  For now, we'll assign 0, 1, -1, 2, -2, 3, -3 ... and so on.
        # We scale this by the point_size*gap_multiplier to get a y coordinate
        # in px.
        row["y"] = (row["y_slot"]//2) * \
            negative_1_if_count_is_odd(row["y_slot"])*point_size*gap_multiplier

    # if the number of points is even, move y-coordinates down to put an equal
    # number of entries above and below the axis
    for row in list_of_rows:
        if bin_counter[row["bin"]] % 2 == 0:
            row["y"] -= point_size*gap_multiplier/2

    df = pd.DataFrame(list_of_rows)

    return df
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


if __name__ == '__main__':

    st.set_page_config(page_title='Dose Dashboard', layout='wide')
    st.title('Patient Dose Dashboard')

    directory = '../../Mined_Data'
    dose_metrics = ['Mean' ,'Median', 'Max', 'Min']

    summary_metrics, number_of_patients = read_csv(directory)

    for metric in dose_metrics:
        summary_metrics[metric] /= 100

    with st.sidebar:

        st.markdown(f'# Patient')
        selected_patient = st.selectbox(label=f'Please select patient to display ({number_of_patients})', options=summary_metrics['Patient'].unique())
        with st.expander('Patient OAR Selection'):
            specific_patient_df = summary_metrics[(summary_metrics['Patient'] == selected_patient)]
            oar_list = specific_patient_df['OAR'].unique().tolist()
            oar_selection = pd.DataFrame({"OAR": oar_list,"Display": [True]*len(oar_list),})
            updated_selection = st.data_editor(oar_selection,
                                                column_config={"Display": st.column_config.CheckboxColumn("Display?",default=True,)},
                                                hide_index=True)
        oar_list = updated_selection['OAR'][updated_selection['Display'] == True]
        st.markdown(f'# Population')
        age_range = st.slider("Pupulation age range", 0, 130, (0,130))

    summary_metrics = summary_metrics[summary_metrics['OAR'].isin(oar_list)]
    
    for oar in oar_list:
        column1, column2 = st.columns([1,1])
        with column1:
            column1_1, column1_2, column1_3 = st.columns([2,1,1])
            with column1_1:
                st.markdown(f'## {oar}')
            with column1_2:
                selected_metric = st.selectbox(label='Please select dose metric', options=dose_metrics, key=oar)
            st.plotly_chart(plot_single_violin_plotly(summary_metrics, specific_patient_df, selected_metric, oar),
                            config={'displayModeBar': False})
            #st.plotly_chart(swarm(summary_metrics[summary_metrics['OAR']==oar][selected_metric], 'test'))
            #st.pyplot(plot_single_violin(summary_metrics, specific_patient_df, selected_metric, oar))


        # with column2:
        #     with st.expander('NTCP Models'):
        #         st.write('Some more')
        #         st.write('Some more')
        #         st.write('Some more')
        #         st.write('Some more')
        #         st.write('Some more')
        #         st.write('Some more')

#st.pyplot(sns.pairplot(summary_metrics.drop('Patient', axis=1), hue='OAR', palette='gnuplot'))
#st.plotly_chart(px.parallel_categories(summary_metrics))

#https://stackoverflow.com/questions/8230638/parallel-coordinates-plot-in-matplotlib


#https://plotly.com/python/parallel-categories-diagram/ (Parallel Categories Linked Brushing)
#https://masamasace.github.io/plotly_color/


