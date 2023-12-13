FROM continuumio/miniconda3

WORKDIR /usr/src/pMTnet_Omni_Document

COPY . .

RUN pip install --no-cache-dir -r requirements.txt

RUN mkdir -p /usr/src/user_input/

