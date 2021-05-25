using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;

namespace InterpolationLab
{
    public partial class Mainform : Form
    {
        public Mainform()
        {
            InitializeComponent();
            Interpolator.DrawGenFunction(chart);
            Interpolator.DrawLagrangeFunction(chart);
            Interpolator.DrawNewtonFunction(chart);
            Interpolator.DrawStandartFunction(chart);
            //Interpolator.DrawSplineFunction(chart);
        }
    }
}
