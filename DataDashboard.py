import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from streamlit_plotly_events import plotly_events
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
import os
import re
import collections


@st.cache_data
def read_database(location):
    path = os.path.join(location, 'All_Patients_Summarised.csv')
    df = pd.read_csv(path)
    return df


@st.cache_data
def get_list_of_metrics_for_contours(location):
    path = os.path.join(location, 'structures_and_metrics.csv')
    df = pd.read_csv(path)
    already_captured = ['max', 'mean', 'median', 'min', 'structurevol']
    df.drop(df[df.Metric.isin(already_captured)].index, inplace=True)
    df['Combined_Metric'] = df.apply(lambda row:f'{row['Metric']}~{row['Value']}', axis=1)
    metrics_for_contours_dict = dict((name,df[df['Structure_Names']==name]['Combined_Metric'].tolist()) for name in df['Structure_Names'].unique().tolist())

    translate_to_human = dict((name,translate(name)) for name in df['Combined_Metric'].unique().tolist() + ['Mean' ,'Median', 'Max', 'Min'])
    translate_from_human = dict((translate(name),name) for name in df['Combined_Metric'].unique().tolist() + ['Mean' ,'Median', 'Max', 'Min'])

    return metrics_for_contours_dict, translate_to_human, translate_from_human


def translate(full_name):
    match full_name:
        case 'Max':
            return 'Max (cGy)'
        case 'Mean':
            return 'Mean (cGy)'
        case 'Median':
            return 'Median (cGy)'
        case 'Min':
            return 'Min (cGy)'

    first = full_name.split('_')[0]
    rest = full_name.split('_')[1]
    second = rest.split('~')[0]
    value = rest.split('~')[1]
    if first == 'dose':
        if second == 'cc':
            formatted_name = f'D{value} cm\u00b3 (cGy)'
        else:
            formatted_name = f'D{value} % (cGy)'
    elif first == 'absolute':
        formatted_name = f'V{value} cGy (cm\u00b3)'
    elif first == 'relative':
        formatted_name = f'V{value} cGy (%)'

    return formatted_name


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


def plot_single_violin_plotly(all_data, patient, patient_id, metric, oar):
    specific_oar_df = all_data[all_data['OAR']==oar]
    violin = go.Violin(x=specific_oar_df[metric],
                        points=False,
                        # jitter=1,
                        # pointpos=0,
                        marker=dict(color='#4575b4'),
                        spanmode='hard',
                        bandwidth=50,
                        name='violin',
                        hoverinfo='none',
                        zorder=2)
    scatter = go.Scatter(x=patient[patient['OAR']==oar][metric],
                            y=['violin'],
                            marker=dict(color='#a50026'),
                            text=patient_id,
                            hovertemplate=
                            '<b>PID</b>: '+ str(patient_id) + '<br>' +
                            '<b>Dose (Gy)</b>: %{x}',
                            hoverlabel = dict(font=dict(color='#a50026')),
                            name='Selected Patient',
                            zorder=3)
    points_df = swarm_df(all_data[all_data['OAR']==oar][['Patient', metric]], metric)
    number_of_patients_for_this_contour = len(points_df)
    points = go.Scatter(x=points_df['x'],
                            y=points_df['y'],
                            marker=dict(color='#abd9e9', symbol='circle-open'),
                            mode='markers',
                            yaxis='y2',
                            text=points_df['patient'],
                            hovertemplate=
                            '<b>PID</b>: %{text}<br>' +
                            '<b>Dose (Gy)</b>: %{x}',
                            name='All Patients',
                            zorder=1)
    
    highlight_df = swarm_df(all_data[(all_data['OAR']==oar)][['Patient', metric]], metric)
    highlight_df = highlight_df[highlight_df['patient'].isin(st.session_state.selected_patient_ids)]
    highlight = go.Scatter(x=highlight_df['x'],
                            y=highlight_df['y'],
                            marker=dict(color='#f46d43', symbol='circle'),
                            mode='markers',
                            yaxis='y2',
                            text=highlight_df['patient'],
                            hovertemplate=
                            '<b>PID</b>: %{text}<br>' +
                            '<b>Dose (Gy)</b>: %{x}',
                            name='Highlighted Patients',
                            zorder=4)
    fig = go.Figure(data=[highlight,scatter,points,violin])
    
    fig.update_layout(yaxis=dict(visible=False),
                        yaxis2=dict(visible=False, overlaying='y', side='right'),
                        xaxis=dict(showgrid=True),
                        margin=dict(l=20, r=20, t=20, b=20),
                        height=150,
                        showlegend=False)
    return fig, number_of_patients_for_this_contour


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

def swarm_df(
    X_series,
    metric,
    point_size=0.5,
    fig_width=800,
    gap_multiplier=1.1,
    bin_fraction=0.95,  # slightly undersizes the bins to avoid collisions
):
    # sorting will align columns in attractive c-shaped arcs rather than having
    # columns that vary unpredictably in the x-dimension.
    # We also exploit the fact that sorting means we see bins sequentially when
    # we add collision prevention offsets.
    X_series = X_series.copy().sort_values(by=metric)

    # we need to reason in terms of the marker size that is measured in px
    # so we need to think about each x-coordinate as being a fraction of the way from the
    # minimum X value to the maximum X value
    min_x = min(X_series[metric])
    max_x = max(X_series[metric])

    list_of_rows = []
    # we will count the number of points in each "bin" / vertical strip of the graph
    # to be able to assign a y-coordinate that avoids overlapping
    bin_counter = collections.Counter()

    for pid, x_val in X_series.values:
        # assign this x_value to bin number
        # each bin is a vertical strip slightly narrower than one marker
        bin = (((fig_width*bin_fraction*(x_val-min_x))/(max_x-min_x)) // point_size)

        # update the count of dots in that strip
        bin_counter.update([bin])

        # remember the "y-slot" which tells us the number of points in this bin and is sufficient to compute the y coordinate unless there's a collision with the point to its left
        list_of_rows.append(
            {"x": x_val, "y_slot": bin_counter[bin], "bin": bin, 'patient': pid})

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


def set_new_patient(value):
    st.session_state.new_patient = value


def store_selected_ids(selection):
    # if len(st.session_state.selected_patient_ids) == 0:
    # if len(selection['selection']['points']) != 0:
    # if selected_patients['selection']['points'][0]['curve_number'] == 0:
    #     selected_patients['selection']['points'].pop(0)
    # patient_ids = [x['text'] for x in selected_patients['selection']['points']]
    # st.session_state.selected_patient_ids = patient_ids
    patient_ids = []
    for point in selected_patients['selection']['points']:
        if point['y'] != 'violin':
            patient_ids.append(point['text'])
    if st.session_state.selected_patient_ids != patient_ids:
        st.session_state.selected_patient_ids = patient_ids
        st.rerun()


if __name__ == '__main__':

    st.set_page_config(page_title='Dose Dashboard', layout='wide')
    st.title('Patient Dose Dashboard')

    # Initiate new session states
    if 'new_patient' not in st.session_state:
        st.session_state['new_patient'] = True
    if 'selected_patient_ids' not in st.session_state:
        st.session_state['selected_patient_ids'] = []

    directory = './'  #'../../Mined_Data'
    summary_metrics = read_database(directory)
    contour_metrics_dict, tth, tfh = get_list_of_metrics_for_contours(directory)

    number_of_patients = len(summary_metrics['Patient'].unique())

    with st.sidebar:
        st.markdown(f'# Patient')
        selected_patient = st.selectbox(label=f'Please select patient to display ({number_of_patients})',
                                        options=summary_metrics['Patient'].unique(),
                                        on_change=set_new_patient,args=(True,))
        with st.expander('Patient Contours Selection'):
            treatment_planning_contours = ['CTV_High', 'PTV_High']
            button1_loc, button2_loc, button3_loc = st.columns([1,1,1])
            specific_patient_df = summary_metrics[(summary_metrics['Patient'] == selected_patient)]
            oar_list = specific_patient_df['OAR'].unique().tolist()

            if st.session_state.new_patient:
                st.session_state.oar_selection = pd.DataFrame({"OAR": oar_list,"Display": [True]*len(oar_list),})
                st.session_state.new_patient = False
            with button1_loc:
                if st.button('All', width='stretch'):
                    st.session_state.oar_selection['Display'].values[:] = True
            with button2_loc:
                if st.button('Planning', width='stretch'):
                    for name, val in st.session_state.oar_selection.values:
                        if (name in treatment_planning_contours):
                            st.session_state.oar_selection.loc[st.session_state.oar_selection['OAR']==name, 'Display'] = True
                        else:
                            st.session_state.oar_selection.loc[st.session_state.oar_selection['OAR']==name, 'Display'] = False
            with button3_loc:
                if st.button('OARs', width='stretch'):
                    for name, val in st.session_state.oar_selection.values:
                        if (name not in treatment_planning_contours):
                            st.session_state.oar_selection.loc[st.session_state.oar_selection['OAR']==name, 'Display'] = True
                        else:
                            st.session_state.oar_selection.loc[st.session_state.oar_selection['OAR']==name, 'Display'] = False

            updated_selection = st.data_editor(st.session_state.oar_selection,
                                                column_config={"Display": st.column_config.CheckboxColumn("Display?")},
                                                hide_index=True)
        oar_list = updated_selection['OAR'][updated_selection['Display'] == True]
        st.markdown(f'# Population')
        if st.button('Clear Selection'):
            st.session_state.selected_patient_ids = []
        with st.expander('WIP'):
            age_range = st.slider("Pupulation age range", 0, 130, (0,130))
        st.markdown(f'# Models')
        with st.expander('WIP'):
            st.write('Some selection for which NTCP models to display')

    summary_metrics = summary_metrics[summary_metrics['OAR'].isin(oar_list)]
    
    for oar in oar_list:
        try:
            contour_specific_dose_metrics = contour_metrics_dict[oar]
        except:
            contour_specific_dose_metrics = []
        dose_metrics = ['Mean' ,'Median', 'Max', 'Min'] + contour_specific_dose_metrics
        dose_metrics_to_display = [tth[x] for x in dose_metrics]
        column1, column2 = st.columns([1,1])
        with column1:
            column1_1, column1_2 = st.columns([2,1])
            with column1_2:
                selected_metric = st.selectbox(label='Please select dose metric', options=dose_metrics_to_display, key=oar)
                selected_metric = tfh[selected_metric]
                fig, this_contour_patient_number = plot_single_violin_plotly(summary_metrics, specific_patient_df, selected_patient, selected_metric, oar)
            with column1_1:
                st.markdown(f'## {oar} ({this_contour_patient_number})')
            selected_patients = st.plotly_chart(fig, on_select='rerun', config={'modeBarButtonsToAdd':['select2d','lasso2d'],'modeBarButtonsToRemove': ['pan',  'zoomIn', 'zoomOut']})
            if len(selected_patients['selection']['points']) != 0:
                store_selected_ids(selected_patients)
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


