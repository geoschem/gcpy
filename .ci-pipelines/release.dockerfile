FROM conda/miniconda3

COPY . /gcpy
RUN conda update -n base -c defaults conda
RUN conda install --file /gcpy/requirements.txt -c conda-forge
RUN apt-get update && apt-get install -y vim
RUN echo "#!/usr/bin/env bash" > /usr/bin/start-container.sh \
&&  echo 'export PYTHONPATH=/gcpy' >> /usr/bin/start-container.sh \
&&  echo 'if [ $# -gt 0 ]; then exec "$@"; else /bin/bash ; fi' >> /usr/bin/start-container.sh \
&&  chmod +x /usr/bin/start-container.sh
ENTRYPOINT ["start-container.sh"]
