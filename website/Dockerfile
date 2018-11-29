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

COPY requirements.txt /
RUN virtualenv -p /usr/bin/python3 /venv/; . /venv/bin/activate && pip3 install -r requirements.txt;
RUN . /venv/bin/activate && python -m ipykernel install --user --name venv --display-name "python3-venv";

COPY boot.sh wsgi.py /website/
COPY app /website/app
RUN chmod +x /website/boot.sh; rm -r /website/app/static/data; ln -s /data /website/app/static/data; chmod -R 777 /root;

ENTRYPOINT cd website; ./boot.sh