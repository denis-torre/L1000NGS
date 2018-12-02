FROM ubuntu:17.10

RUN apt-get update && apt-get install -y python3 \
 python3-pip \
 python3-dev \
 python3-setuptools \
 nginx \
 uwsgi-core \
 libpcre3 \
 libpcre3-dev \
 virtualenv \
 uwsgi-plugin-python3

RUN mkdir -p /website/app /data
WORKDIR /website

COPY requirements.txt .
RUN virtualenv -p /usr/bin/python3 venv/; . venv/bin/activate && pip3 install -r requirements.txt;
RUN . venv/bin/activate && python -m ipykernel install --user --name venv --display-name "python3-venv";

COPY data /data
COPY boot.sh wsgi.py ./
COPY app ./app
RUN chmod +x boot.sh; rm -r app/static/data; ln -s /data app/static/data;

ENTRYPOINT ./boot.sh