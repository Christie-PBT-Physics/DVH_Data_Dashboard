import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
import os
import re


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


def plot_single_violin_plotly(all, patient, metric, oar):
    specific_oar_df = all[all['OAR']==oar]
    violin = go.Violin(x=specific_oar_df[metric],
                        points='all',
                        jitter=1,
                        pointpos=0,
                        marker=dict(size=8, symbol="x-thin-open"),
                        spanmode='hard',
                        name='violin')
    scatter = go.Scatter(x=patient[patient['OAR']==oar][metric],
                            y=['violin'],
                            marker=dict(color='red'))
    fig = go.Figure(data=[violin,scatter])
    
    fig.update_layout(yaxis=dict(visible=False),
                        xaxis=dict(showgrid=True),
                        margin=dict(l=20, r=20, t=20, b=20),
                        height=150,
                        showlegend=False)
    return fig
    

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

#https://stackoverflow.com/questions/8230638/parallel-coordinates-plot-in-matplotlib

