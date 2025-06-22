# ARG-Classifier-Nucleotide-Transformer
This repository contains the code and model files for a deep learning-based approach to classify antibiotic resistance genes (ARGs) using a transformer model adapted for nucleotide sequences. The goal of this project is to enable efficient and accurate prediction of ARGs from raw DNA sequences.

Authors
	•	Yesasvi Sai Nandigam
	•	Prasanth Kumar Thuthika

Project Description

Antimicrobial resistance is a growing global concern, and early identification of resistance genes is crucial for both research and clinical decision-making. This project fine-tunes a transformer model previously trained on nucleotide data for binary classification tasks involving ARG detection.

By training on labeled genomic datasets—comprising resistant genes, non-resistant genes, and synthetic sequences—the model learns to distinguish ARGs in diverse sequence contexts. This tool can assist in genomic surveillance and resistome profiling efforts by providing rapid ARG predictions from sequence data alone.

Installation

To run this project, begin by setting up a virtual environment:
conda create -n arg_classifier_env python=3.10
conda activate arg_classifier_env

Install the necessary Python packages:
pip install transformers==4.38.1 datasets==2.18.0 torch==2.2.1 pandas==2.2.1 \
            numpy==1.26.4 scikit-learn==1.4.1 matplotlib==3.8.3 seaborn==0.13.2

Downloading the Model

The fine-tuned transformer model used in this project is available for download from the following link:

ARG Classifier - Google Drive

After downloading, extract the contents and place the arg_classifier/ directory in the root folder of this repository. The model directory should include the following files:
	•	config.json
	•	model.safetensors
	•	tokenizer_config.json
	•	special_tokens_map.json
	•	vocab.txt

Running Predictions

Use the following script to classify a given DNA sequence:

import torch
from transformers import AutoTokenizer, AutoModelForSequenceClassification

Load model and tokenizer
model_directory = "arg_classifier"
tokenizer = AutoTokenizer.from_pretrained(model_directory)
model = AutoModelForSequenceClassification.from_pretrained(model_directory)

Input sequence
sequence = "ATGCGTACGTAGCTAGCTAGCTAGCGTATCGTAGCTAGT"
encoded = tokenizer(sequence, return_tensors="pt", padding="max_length", truncation=True, max_length=512)

Generate prediction
model.eval()
with torch.no_grad():
    outputs = model(**encoded)
    label = torch.argmax(outputs.logits, dim=1).item()

Print result
print("Prediction:", "Resistant" if label == 1 else "Non-Resistant")

Example output:
Prediction: Resistant
