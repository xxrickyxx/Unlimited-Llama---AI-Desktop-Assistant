# 🦙 Unlimited Llama - AI Desktop Assistant AiLo Core

<div align="center">

<p align="center">
  <img src="icon.png" alt="Unlimited Llama" width="128">
</p>

**A complete AI desktop assistant with chat, web search, speech synthesis, and OCR.**

[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Windows](https://img.shields.io/badge/Windows-Supported-success)](https://www.microsoft.com/windows)
[![Open Source](https://badges.frapsoft.com/os/v1/open-source.svg?v=103)](https://opensource.org/)
[![Download on Hugging Face](http://huggingface.co/front/assets/huggingface_logo-noborder.svg)](https://huggingface.co/xxrickyxx/Unlimited-Llama/blob/main/README.md)

</div>

---

## ✨ Features

- 💬 **Smart chat** with local GGUF models  
- 🌐 **Integrated web search** for up-to-date information  
- 🔊 **Text-to-Speech (TTS)** and **Speech Recognition (STT)**  
- 📷 **OCR** to extract text from images  
- 💾 **Advanced session management**  
- 🎛️ **Supports any LLM model size**  
- 🔌 **OpenAI-compatible API server**  
- 📤 **Export** in JSON, TXT, and Markdown  

---

## 📦 Download

[![Download on Hugging Face](https://img.shields.io/badge/Download-HuggingFace-yellow?logo=huggingface)](https://huggingface.co/xxrickyxx/Unlimited-Llama/blob/main/README.md)

---

## 🚀 Quick Start Guide

**First Launch**
1. Load a model → 🤖 *Model → 📁 Load Model*  
2. Start chatting → type in the box below and press *Enter*  
3. Sessions are saved automatically  

---

### 🔍 Web Search
- Enable/disable using the 🌐 *Web Search* toggle  
- Automatically searches for news, recent info, or local data  
- Displays the sources used  

### 🔊 Speech Synthesis (TTS)
- Enable via 🔊 *TTS* in the sidebar  
- The assistant reads responses aloud  
- Use 🔇 *STOP* to interrupt  

### 🎤 Speech Recognition
- 🎤 *Voice Input* for single input  
- 🎤 *Start Listening* for continuous mode  

### 📷 OCR from Images
- Click 📷 *Image OCR*  
- Select an image (PNG, JPG, etc.)  
- Extracted text is automatically inserted into the chat  

---

## 🛠️ Troubleshooting

**❌ “Model not found”**  
- Make sure the GGUF file is in the `/models` folder  
- Verify the file format is `.gguf`  
- Check that you have enough disk space  

**❌ “Tesseract not found”**  
- Install Tesseract OCR following the instructions below  
- Restart the application after installation  

---
⚙️ Configuration
Memory Optimization
Memory Mapping (MMAP)
What it does: Maps model directly from disk instead of loading entirely into RAM

Benefits: Reduces RAM usage by up to 70%, faster startup

Use when: Limited RAM, large models (>7GB)

Performance: Slightly slower inference, much less RAM usage

Memory Locking (MLOCK)
What it does: Locks model in RAM preventing swap to disk

Benefits: Maximum performance, consistent response times

Use when: Abundant RAM, performance-critical applications

Performance: Fastest inference, permanent RAM occupation
---

---

## 🔢 Goldbach Conjecture – Numerical Experiments

This repository includes a standalone numerical experiment tool for the
**fiber-sum / discrete Fourier-shift** approach to the Goldbach conjecture.

### Quick start

```bash
pip install numpy
python goldbach_experiment.py --N 100000 --W 256 --a 1 --q 6 --beta 0.0
```

### Random-beta sweep (50 samples)

```bash
python goldbach_experiment.py --N 100000 --W 256 --a 1 --q 6 --samples 50
```

### Self-check

```bash
python goldbach_experiment.py --selfcheck
```

### Run unit tests

```bash
pip install pytest
python -m pytest tests/test_goldbach_experiment.py -v
```

See [`docs/goldbach_experiment.md`](docs/goldbach_experiment.md) for the full
experimental protocol, mathematical background, and CLI reference.

---

## ⚙️ System Requirements

### Minimum
- **OS:** Windows 10/11, macOS 10.15+, Linux (Ubuntu 18.04+)  
- **RAM:** 8 GB (16 GB recommended)  
- **Disk Space:** 2 GB + space for models  
- **CPU:** Modern 64-bit processor  

### Recommended
- **RAM:** 16 GB+ for large models  
- **GPU:** NVIDIA/AMD with CUDA or Metal (optional)  
- **Disk Space:** 10 GB+ for large models  

---

## 🔧 Installation

### 1. Install Tesseract OCR (Required for OCR)

#### Windows
```bash
# Using Chocolatey (recommended)
choco install tesseract
