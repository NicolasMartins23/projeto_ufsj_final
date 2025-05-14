import streamlit as st

from BioApps import DNA
from extras import SEQUENCIA_DNA, Downloads

# Coloque seu nome aqui
name = 'Nicolas Martins'

st.header("🖥️🧬 Meu 1º BioApp 🧬🖥️")
st.write("Aplicação para análise de Bioinformática de sequências de DNA")

# input
dna_input = st.text_input(label="Inserir aqui a sequência de DNA",
                          value=SEQUENCIA_DNA)

button_input = st.button("Analisar")

if button_input:
    my_dna = DNA(dna_input)
    st.subheader("Sequência de DNA:")
    st.code(my_dna.sequence)
    my_rna = my_dna.to_rna()
    st.subheader("Sequência de RNA:")
    st.code(my_rna.sequence)
    my_protein = my_rna.to_protein()
    st.subheader("Sequência Peptídica:")
    st.code(my_protein.sequence)

    protein_info = [
        {
            "Propriedade": "Sequência Peptídica",
            "Valor": my_protein.sequence
        },
        {
            "Propriedade": "Massa Molecular",
            "Valor": my_protein.molecular_weight(),
        },
    ]

    st.subheader("Resultados")
    st.table(protein_info)

    # Cria o objeto de gerenciamento de downloads
    downloads = Downloads()

    file_name_input = st.text_input(label="Nome do Arquivo",
                                    value="Results_01")
    protein_name_input = st.text_input(label="Nome da Proteína",
                                       value="Protein_01")

    downloads.zip_download(file_name=file_name_input,
                           csv_results=protein_info,
                           protein_name=protein_name_input,
                           protein_sequence=my_protein.sequence)

st.write("---")
st.write(f"Desenvolvido por {name}")
