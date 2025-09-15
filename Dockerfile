FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV LIBGL_ALWAYS_INDIRECT=1

# Add user
RUN adduser --quiet --disabled-password qtuser

# System dependencies for Python + Qt GUI
RUN apt-get update && apt-get install -y \
    software-properties-common \
    python3 \
    python3-pip \
    python3-distutils \
    libgl1 \
    libegl1 \
    libglu1-mesa \
    mesa-utils \
    libxkbcommon-dev \
    # Qt xcb dependencies  \
    libx11-xcb1 \
    libxcb1 \
    libxcb-xinerama0 \
    libxcb-cursor0 \
    libxcb-randr0 \
    libxcb-shape0 \
    libxcb-xfixes0 \
    libxcb-shm0 \
    libxcb-sync1 \
    libxcb-util1 \
    libxcb-keysyms1 \
    libxcb-icccm4 \
    libxcb-render-util0 \
    libxcb-image0 \
    libfontconfig1 \
    libfreetype6 \
    libxrender1 \
    libxext6 \
    libx11-6 \
    libxi6 \
    libxcomposite1 \
    libxcursor1 \
    libxdamage1 \
    libxrandr2 \
    libxtst6 \
    # Ensure full Qt platform plugin deps \
    qt6-base-dev \
    qt6-base-dev-tools \
    && rm -rf /var/lib/apt/lists/*

# Point Qt to plugin directory
ENV QT_QPA_PLATFORM_PLUGIN_PATH=/usr/lib/x86_64-linux-gnu/qt6/plugins/platforms

# Create working directories
RUN mkdir /HOST && mkdir -p /home/software/Icons
WORKDIR /home/software

# Python dependencies
RUN pip install --no-cache-dir \
    cellpose==3.1.0 \
    PyQt6==6.8.1 \
    pyqtgraph==0.13.7 \
    tifffile==2023.2.28 \
    XlsxWriter==3.2.2 \
    scikit-learn==1.6.1 \
    scikit-image==0.24.0 \
    aicsimageio==4.14.0

# Copy project files
COPY keys_size_factor.npy mycmap.bin /home/software/
COPY *.py /home/software/
COPY Icons/* /home/software/Icons/

# Permissions
RUN chmod -R 777 /home/software

# Default command
CMD ["python3", "TrojanHorseQuantifier_v1_0.py"]




# INFO HOW TO RUN
# https://hub.docker.com/r/fadawar/docker-pyqt5
# docker build -t first_gui .
# docker run --rm -it   -v /home/atrullo:/HOST  -v /tmp/.X11-unix:/tmp/.X11-unix  -e DISPLAY=$DISPLAY   -u qtuser   trojan_horse_quantifier  #  THIS WORKS ON LOCAL MACHINE
# docker run --rm -it   -v /home/atrullo:/HOST  -v /tmp/.X11-unix:/tmp/.X11-unix  -e DISPLAY=$DISPLAY   follow_cells   #  THIS WORKS ON LOCAL MACHINE
# docker run --net=host --env="DISPLAY" -v /home/fabrice:/HOST --volume="$HOME/.Xauthority:/root/.Xauthority:rw" track_split    #  THIS WORKS ON THE SERVER INSTEAD (DOCKER ON VIRTUAL MACHINE ON A SERVER)
