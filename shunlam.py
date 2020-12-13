import numpy as np
from scipy.io.wavfile import write
from pydub import AudioSegment
from pydub.playback import play


data = np.random.uniform(-1,1,44100) # 44100 random samples between -1 and 1
scaled = np.int16(data/np.max(np.abs(data)) * 32767)
write('sound.wav', 44100, scaled)
song = AudioSegment.from_wav("sound.wav")
play(song)
