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

## 🔬 Experimental: Goldbach Fiber-Sum Tools (S1.2 / Option 2)

The `tools/` directory contains numerical experiments for the Goldbach
fiber-sum approach.  Full details are in
[`docs/goldbach/experiment-plan.md`](docs/goldbach/experiment-plan.md).

### Requirements
```bash
pip install numpy   # optional but recommended for speed
```

### CLI Usage

```bash
# Self-check / doctests (no dependencies required)
python tools/goldbach_fiber_experiments.py --selfcheck

# Basic run: N=10000, W=30, base α=0/1 (major arc), 10 minor-arc samples
python tools/goldbach_fiber_experiments.py --N 10000 --W 30

# Minor-arc experiment: α near 3/7, W=210, 200 samples
python tools/goldbach_fiber_experiments.py --N 50000 --W 210 --a 3 --q 7 --samples 200

# Convergence sweep: N doubles from 1000 to 32000, W=30
python tools/goldbach_fiber_experiments.py --N 32000 --W 30 --sweep --N0 1000

# Full experiment with custom major-arc threshold Q
python tools/goldbach_fiber_experiments.py --N 100000 --W 2310 --a 1 --q 11 --samples 500 --Q 0.01

# Quiet mode (suppress per-sample output, print summary only)
python tools/goldbach_fiber_experiments.py --N 20000 --W 30 --samples 20 --quiet
```

### Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--N` | 10000 | Prime upper bound |
| `--W` | 30 | Modulus (try primorials: 6, 30, 210, 2310) |
| `--a` | 0 | Numerator of base α = a/q |
| `--q` | 1 | Denominator of base α = a/q |
| `--samples` | 10 | Number of minor-arc α values to probe |
| `--Q` | auto | Major-arc half-width threshold |
| `--sweep` | off | Enable N-convergence sweep (M5) |
| `--N0` | 1000 | Starting N for sweep |
| `--selfcheck` | off | Run doctests and exit |
| `--quiet` | off | Suppress verbose per-sample output |

---

## 🔧 Installation

### 1. Install Tesseract OCR (Required for OCR)

#### Windows
```bash
# Using Chocolatey (recommended)
choco install tesseract
