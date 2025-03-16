import numpy as np
import pandas as pd
import streamlit as st
import joblib  # Use joblib instead of pickle
from sklearn.preprocessing import StandardScaler

# Load the model with error handling
try:
    load_model = joblib.load('Final_model.joblib')
    st.success("Model loaded successfully!")
except Exception as e:
    st.error(f"Error loading model: {e}")
    st.stop()  # Stop the app if the model fails to load

# Function to process input sequence
def input_seq():
    if 'sequence1' not in st.session_state:
        st.session_state.sequence1 = ""

    str22 = "CGCCUCCCACGCGGGAGACCCGGGUUCAAUUCCCGGCCAAU"

    if st.button('Sample Sequence'):
        st.session_state.sequence1 = str22

    sequence1 = st.text_area("Sequence Input", value=st.session_state.sequence1, height=200)
    st.session_state.sequence1 = sequence1

    if st.button("Submit"):
        abc = sequence1
        count = []
        keeper = []
        len1 = len(abc)
        len2 = len1 - 1

        for i in range(len1):
            if abc[i] == "C":
                count.append(i)

        len3 = len(count)
        for i in range(len3):
            s = count[i]
            if s <= 20 and s <= len2:
                n = len2 - s
                m = 20 - s
                str1 = ("U" * m)
                str2 = abc[s - s:s]

                if n <= 20:
                    o = 20 - n
                    str4 = abc[s:s + n + 1]
                    str5 = ("U" * o)
                    str6 = "".join((str1, str2, str4, str5))
                    keeper.append(str6)
                elif n > 20:
                    str4A = abc[s:s + 21]
                    str6A = "".join((str1, str2, str4A))
                    keeper.append(str6A)

            elif s > 20:
                n1 = len1 - s
                str7 = abc[s - 20:s + 1]
                if n1 <= 20:
                    o1 = 20 - n1
                    o1 = o1 + 1
                    str9 = ("C" * o1)
                    str8 = abc[s:s + n1 - 1]
                    str10 = "".join((str7, str8, str9))
                    keeper.append(str10)
                elif n1 > 20:
                    str8A = abc[s:s + 20]
                    str10A = "".join((str7, str8A))
                    keeper.append(str10A)

        klen = len(keeper)
        allFVs = []
        for i in range(klen):
            seq = keeper[i]
            allFVs.append(calcFV(seq.lower()))

        with open('./IISequence_FVs_for_test.csv', mode='w') as fvFile:
            fvWriter = csv.writer(fvFile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            for fv in allFVs:
                fvWriter.writerow(fv)

        df = pd.read_csv("IISequence_FVs_for_test.csv", sep=',', header=None)
        W = df.iloc[:, :].values
        std_scale = StandardScaler().fit(W)
        W = std_scale.transform(W)

        pred = load_model.predict(W)
        output_proba = load_model.predict_proba(W)[:, 1]

        lno = len(output_proba)
        for i in range(lno):
            st.subheader("Site Number = ")
            st.write(count[i])
            st.subheader("Sequence")
            st.write(keeper[i])

            if output_proba[i] > 0.7:
                st.info("Output = 5-Methylcytosine Site")
            else:
                st.info("Output = Non-5-Methylcytosine Site")

# Run the function
input_seq()
